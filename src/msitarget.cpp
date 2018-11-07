#include "msitarget.h"
#include <iostream>
#include <fstream>
#include "util.h"
#include <string.h>
#include <sstream>
#include "bamutil.h"
#include <math.h>
#include "emd.h"

MsiTarget::MsiTarget(string chr, int start, int end, string name){
    mChr = chr;
    mStart = start;
    mEnd = end;
    mCenter = (mStart + mEnd) /2;
    mName = name;
    mHistogram = new int[MAX_HIST_BIN];
    mMaxBin = 0;
    mMinBin = MAX_HIST_BIN-1;
    mSupportingReads = 0;
    memset(mHistogram, 0, sizeof(int)*MAX_HIST_BIN);
}

MsiTarget::~MsiTarget() {
    delete[] mHistogram;
    map<string, Pair*>::iterator iter;
    for(iter = mPairs.begin(); iter!=mPairs.end(); iter++)
        delete iter->second;
}

void MsiTarget::addLen(int len) {
    if(len >MAX_HIST_BIN)
        return;
    mHistogram[len]++;
    mSupportingReads++;
    if(len > mMaxBin)
        mMaxBin = len;
    if(len < mMinBin)
        mMinBin = len;
}

void MsiTarget::addRead(bam1_t *b) {
    string qname = BamUtil::getQName(b);
    if(mPairs.count(qname))
        mPairs[qname]->setRight(b);
    else {
        mPairs[qname] = new Pair();
        mPairs[qname]->setLeft(b);
    }
}

double MsiTarget::emdWith(MsiTarget* other) {
    if(mSupportingReads == 0 || other->mSupportingReads==0)
        return 0.0;
    int minBin = min(mMinBin, other->mMinBin);
    int maxBin = max(mMaxBin, other->mMaxBin);
    int bins = maxBin - minBin + 1;
    double* thisWeights = new double[bins];
    double* otherWeights = new double[bins];
    for(int i=0; i<bins; i++) {
        thisWeights[i] = mHistogram[i+minBin] / (double)mSupportingReads;
        otherWeights[i] = other->mHistogram[i+minBin] / (double)other->mSupportingReads;
    }
    double** costs = new double*[bins];
    double** flows = new double*[bins];
    for(int i=0; i<bins; i++) {
        costs[i] = new double[bins];
        flows[i] = new double[bins];
    }
    for(int i=0; i<bins; i++) {
        for(int j=0; j<bins; j++) {
            costs[i][j] = abs(i-j);
        }
    }

    double emdValue = emd(bins, thisWeights, bins, otherWeights, costs, flows);

    for(int i=0; i<bins; i++) {
        delete[] costs[i];
        delete[] flows[i];
    }

    delete[] flows;
    delete[] costs;
    delete[] thisWeights;
    delete[] otherWeights;
    return emdValue;
}

double MsiTarget::entropy() {
    if(mSupportingReads == 0)
        return 0;
    double total = (double)mSupportingReads;

    double entropyValue = 0.0;
    for(int i=mMinBin; i<=mMaxBin; i++) {
        if(mHistogram[i] > 0) {
            double p = (mHistogram[i] / total);
            entropyValue += -p*log2(p);
        }
    }
    return entropyValue;
}

vector<MsiTarget*> MsiTarget::parseBed(string filename) {
    ifstream file;
    file.open(filename.c_str(), ifstream::in);
    const int maxLine = 4096;
    char line[maxLine];
    vector<MsiTarget*> targets;
    while(file.getline(line, maxLine)){
        // trim \n, \r or \r\n in the tail
        int readed = strlen(line);
        if(readed >=2 ){
            if(line[readed-1] == '\n' || line[readed-1] == '\r'){
                line[readed-1] = '\0';
                if(line[readed-2] == '\r')
                    line[readed-2] = '\0';
            }
        }
        string linestr(line);
        linestr = trim(linestr);
        vector<string> splitted;
        split(linestr, splitted, "\t");
        // wrong line
        if(splitted.size()<2)
            continue;
        // comment line
        if(starts_with(splitted[0], "#"))
            continue;
        // position line require id, start, position
        if(splitted.size()<4)
            continue;

        string chr = trim(splitted[0]);
        int start = atoi(trim(splitted[1]).c_str());
        int end = atoi(trim(splitted[2]).c_str());
        string name = trim(splitted[3]);
        MsiTarget* t = new MsiTarget(chr, start, end, name);
        targets.push_back(t);
    }
    return targets;
}

void MsiTarget::stat() {
    map<string, Pair*>::iterator iter;
    for(iter = mPairs.begin(); iter!=mPairs.end(); iter++) {
        string inserted = iter->second->getInsertedByPCR(mLeftAdapter, mRightAdapter);
        if(inserted.length()>0) {
            addLen(inserted.length());
        }
        /*
        cerr << mLeftAdapter<<"\t"<<mRightAdapter<<endl;
        if(iter->second->mLeft)
            cerr << BamUtil::getSeq(iter->second->mLeft) << endl;
        if(iter->second->mRight)
            cerr << BamUtil::getSeq(iter->second->mRight) << endl;
        cerr << "inserted: " << inserted << endl;
        cerr << endl;*/
    }
}

void MsiTarget::print() {
    cerr << mChr << "\t" << mStart << "\t" << mEnd << "\t" << mName << "\t" << mLeftAdapter << "\t" << mInserted << "\t" << mRightAdapter << endl;
}
