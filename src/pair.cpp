#include "pair.h"
#include "bamutil.h"
#include <memory.h>
#include "editdistance.h"

Pair::Pair(){
    mLeft = NULL;
    mRight = NULL;
}

Pair::~Pair(){
    if(mLeft) {
        bam_destroy1(mLeft);
        mLeft = NULL;
    }
    if(mRight) {
        bam_destroy1(mRight);
        mRight = NULL;
    }
}

void Pair::setLeft(bam1_t *b) {
    mLeft = b;
    mLeftCigar = BamUtil::getCigar(mLeft);
}

void Pair::setRight(bam1_t *b) {
    mRight = b;
    mRightCigar = BamUtil::getCigar(mRight);
}

bool Pair::pairFound() {
    return mLeft != NULL && mRight != NULL;
}

int Pair::getLeftRef() {
    if(mLeft == NULL)
        return -1;

    return mLeft->core.tid;
}

int Pair::getLeftPos() {
    if(mLeft == NULL)
        return -1;

    return mLeft->core.pos;
}

int Pair::getRightRef() {
    if(mRight == NULL)
        return -1;

    return mRight->core.tid;
}

int Pair::getRightPos() {
    if(mRight == NULL)
        return -1;

    return mRight->core.pos;
}

int Pair::getTLEN() {
    if(mLeft != NULL)
        return abs(mLeft->core.isize);
    else if(mRight != NULL)
        return abs(mRight->core.isize);
    else
        return -1;
}

string Pair::getQName() {
    if(mLeft != NULL)
        return BamUtil::getQName(mLeft);
    else if(mRight != NULL)
        return BamUtil::getQName(mRight);
    else
        return "";
}

Pair::MapType Pair::getMapType() {
    if(mLeft == NULL || mRight == NULL)
        return Unknown;

    int lref = getLeftRef();
    int rref = getRightRef();

    if(lref == rref) {
        if(lref >= 0)
            return ProperlyMapped;
        else
            return NoneMapped;
    }

    if(lref<0 && rref>=0)
        return OnlyRightMapped;

    if(lref>=0 && rref<0)
        return OnlyLeftMapped;

    if(lref != rref)
        return CrossRefMapped;

    return Unknown;
}

string Pair::getLeftCigar() {
    return mLeftCigar;
}

string Pair::getRightCigar() {
    return mRightCigar;
}

string Pair::getInsertedByPCR(string leftAdapter, string rightAdapter) {
    if(!mLeft && !mRight)
        return string();

    // single-end
    if(!mLeft)
        return getInsertedByPCR(BamUtil::getSeq(mRight), leftAdapter, rightAdapter);
    else if(!mRight)
        return getInsertedByPCR(BamUtil::getSeq(mLeft), leftAdapter, rightAdapter);

    // paired-end
    string left = getInsertedByPCR(BamUtil::getSeq(mLeft), leftAdapter, rightAdapter);
    string right = getInsertedByPCR(BamUtil::getSeq(mRight), leftAdapter, rightAdapter);

    if(left.length() > 0 && right.length()>0) { // both reads cover adapters
        // pair with many sequencing errors, skip it
        if(abs((int)left.length() - (int)right.length()) > 3)
            return "";
        int diff, gapPos;
        diff = diffWithOneGap(left.c_str(), right.c_str(), min(left.length(), right.length()), gapPos);
        // pair with many sequencing errors, skip it
        if(diff > 8)
            return "";
        else if(left.length() > right.length()) // we simply use the longer read since DEL is more common than INDEL for sequence erros
            return left;
        else
            return right;
    } else if(left.length() > 0) // only left read covers adapters
        return left;
    else if(right.length() > 0) // only right read covers adapters
        return right;

    // if neighter read covers the adapters, we should merge the this pair
    int read1Left = abs(mLeft->core.pos);
    int read1Right = BamUtil::getLastNonclipRef(mLeft) + read1Left;

    int read2Left = abs(mRight->core.pos);
    int read2Right = BamUtil::getLastNonclipRef(mRight) + read2Left;

    // not overlapped, skip it
    if(read2Left > read1Right) {
        //cerr << "not overlapped " << read2Left << ", " << read1Right << endl;
        return "";
    }

    // otherwise, we use the center to merge this pair
    int center = (read1Right + read2Left)/2;
    int lbreak = -1;
    int rbreak = -1;

    for(int i=0; i<mLeft->core.l_qseq; i++) {
        if(BamUtil::getRefOffset(mLeft, i)+read1Left >= center) {
            lbreak = i;
            break;
        }
    }

    for(int i=0; i<mRight->core.l_qseq; i++) {
        if(BamUtil::getRefOffset(mRight, i)+read2Left >= center) {
            rbreak = i;
            break;
        }
    }

    string leftSeq = BamUtil::getSeq(mLeft);
    string rightSeq = BamUtil::getSeq(mRight);

    //cerr << "merging..." << endl;

    int olen = leftSeq.length() - lbreak + rbreak;
    if(olen > leftSeq.length())
        olen = leftSeq.length();
    if(olen > rightSeq.length())
        olen = rightSeq.length();

    string leftOverlap = leftSeq.substr(leftSeq.length()-olen, olen);
    string rightOverlap = rightSeq.substr(0, olen);

    unsigned int ed = edit_distance(leftOverlap, rightOverlap);
    // this is a bad overlap, should be caused by a false positive merging
    if(ed > 20 || ed*3 > olen) {
        //cerr << "-------- skip merging due to edit distance too high: " << ed << endl;
        //cerr << leftOverlap << endl;
        //cerr << rightOverlap << endl;
        return "";
    }

    string mergedSeq = leftSeq.substr(0, lbreak) + rightSeq.substr(rbreak, rightSeq.length() - rbreak);
    //cerr << "merged with edit distance: " << ed << endl;
    //cerr << mergedSeq << endl;
    return getInsertedByPCR(mergedSeq, leftAdapter, rightAdapter);
}

string Pair::getInsertedByPCR(string seq, string leftAdapter, string rightAdapter) {
    int left=-1;
    int right=-1;
    int slen = seq.length();
    int llen = leftAdapter.length();
    int rlen = rightAdapter.length();
    const char* sdata = seq.c_str();
    const char* ldata = leftAdapter.c_str();
    const char* rdata = rightAdapter.c_str();

    if(slen < llen + rlen)
        return "";

    for(int i=0; i<slen - llen; i++) {
        int diff, gapPos;
        diff = diffWithOneGap(sdata+i, ldata, llen, gapPos);
        if(gapPos == 0) { // no gap
            if(diff <= 1) {
                left = i+llen;
                break;
            }
        } else if(gapPos > 0) { // gap at seq
            if(diff == 0) {
                left = i+llen-1;
                break;
            }
        } else if(gapPos < 0) { // gap at adapter
            if(diff == 0) {
                left = i+llen+1;
                break;
            }
        }
    }

    if(left < 0)
        return "";

    for(int i=0; i<slen - rlen; i++) {
        int diff, gapPos;
        diff = diffWithOneGap(sdata+i, rdata, rlen, gapPos);
        if(diff == 0 || (diff <= 1 && gapPos==0)) {
            right = i-1;
            break;
        }
    }

    if(right<left)
        return "";

    return seq.substr(left, right - left + 1);
}

bool Pair::test() {
    string seq = "ATTATCTGAATATTTAAGGTCTGCCTTAACGTGATCCCCATTGCTGAATTTTACCTCCTGACTCCAAAAACTCTTCTCTTCCCTGGGCCCAGTCCTATTTTTTTTTTTTTTGTGAGACAGAGTCTCACTCTGTCACCCAGGTTGGAATGC";
    string leftAdapter = "CCCCATTGCTGA";
    string rightAdapter = "AGGTTGGAATGC";
    string inserted = getInsertedByPCR(seq, leftAdapter, rightAdapter);
    cerr << "inserted DNA is: " << inserted << endl;
    if(inserted != "ATTTTACCTCCTGACTCCAAAAACTCTTCTCTTCCCTGGGCCCAGTCCTATTTTTTTTTTTTTTGTGAGACAGAGTCTCACTCTGTCACCC")
        return false;
    return true;
}

void Pair::dump() {
    if(mLeft)
        BamUtil::dump(mLeft);
    if(mRight)
        BamUtil::dump(mRight);
}

