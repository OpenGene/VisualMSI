#include "visualmsi.h"
#include "bamutil.h"
#include "jsonreporter.h"
#include "htmlreporter.h"
#include "reference.h"

VisualMSI::VisualMSI(Options *opt){
    mOptions = opt;
}

VisualMSI::~VisualMSI(){
    for(int m=0; m<mMsiTargetsTumor.size(); m++)
        delete mMsiTargetsTumor[m];
    for(int m=0; m<mMsiTargetsNormal.size(); m++)
        delete mMsiTargetsNormal[m];
}

void VisualMSI::calcLocusTarget(vector<MsiTarget*>& msiTargets, map<string, map<int, MsiTarget*>>& locusTarget) {
    for(int m=0; m<msiTargets.size(); m++) {
        string chr = msiTargets[m]->mChr;
        int center = msiTargets[m]->mCenter;
        if(locusTarget.count(chr) == 0)
            locusTarget[chr] = map<int, MsiTarget*>();
        locusTarget[chr][center] = msiTargets[m];
    }
}

void VisualMSI::reportText() {
    bool hasNormal = !mOptions->normal.empty() && mMsiTargetsNormal.size() == mMsiTargetsTumor.size();
    for(int m=0; m<mMsiTargetsTumor.size(); m++) {
        cout<<endl;
        cout << m+1 << ", " << mMsiTargetsTumor[m]->mName << ":" << endl;
        MsiTarget* tumor = mMsiTargetsTumor[m];
        MsiTarget* normal = NULL;
        if(hasNormal)
            normal = mMsiTargetsNormal[m];
        if(hasNormal)
            cout << "the earth mover's distance (EMD): " << tumor->emdWith(normal) <<endl;

        cout<<"entropy of tumor data: " << tumor->entropy() <<endl;
        if(hasNormal)
            cout<<"entropy of normal data: " << normal->entropy() <<endl;

        cout<<"supporting reads of tumor data: " << tumor->mSupportingReads <<endl;
        if(hasNormal)
            cout<<"supporting reads of normal data: " << normal->mSupportingReads <<endl;
        if(tumor->mSupportingReads < mOptions->depthReq || (normal && normal->mSupportingReads < mOptions->depthReq))
            cout << "quality control: failed" << endl;
        else
            cout << "quality control: passed" << endl;
    }
}

void VisualMSI::reportJSON() {
    JsonReporter reporter(mOptions);
    reporter.report(mMsiTargetsTumor, mMsiTargetsNormal);
}

void VisualMSI::reportHTML() {
    HtmlReporter reporter(mOptions);
    reporter.report(mMsiTargetsTumor, mMsiTargetsNormal);
}

void VisualMSI::makeAdapters(vector<MsiTarget*>& msiTargets, map<string, map<int, MsiTarget*>>& locusTarget) {
    const int half = mOptions->targetInsertedSize / 2;
    Reference* ref = Reference::instance(mOptions);
    for(int m=0; m<msiTargets.size(); m++) {
        string chr = msiTargets[m]->mChr;
        string chrSeq = ref->getContig(chr);
        if(chrSeq.empty())
            continue;
        int center = msiTargets[m]->mCenter;
        int leftStart = center - half - mOptions->adapterLen;
        int rightStart = center + half;
        msiTargets[m]->mLeftAdapter = chrSeq.substr(leftStart, mOptions->adapterLen);
        msiTargets[m]->mRightAdapter = chrSeq.substr(rightStart, mOptions->adapterLen);
        msiTargets[m]->mInserted = chrSeq.substr(center - half, half*2);
    }
}

void VisualMSI::stat(vector<MsiTarget*>& msiTargets, map<string, map<int, MsiTarget*>>& locusTarget) {
    for(int m=0; m<msiTargets.size(); m++) {
        msiTargets[m]->stat();
    }
}

void VisualMSI::run() {
    mMsiTargetsTumor = MsiTarget::parseBed(mOptions->targetFile);
    calcLocusTarget(mMsiTargetsTumor, mLocusTargetTumor);
    cerr << "parsing " << mOptions->targetFile << endl;
    makeAdapters(mMsiTargetsTumor, mLocusTargetTumor);
    for(int m=0; m<mMsiTargetsTumor.size(); m++) {
        mMsiTargetsTumor[m]->print();
    }
    evaluate(mOptions->input, mMsiTargetsTumor, mLocusTargetTumor);

    if(!mOptions->normal.empty()) {
        mMsiTargetsNormal = MsiTarget::parseBed(mOptions->targetFile);
        calcLocusTarget(mMsiTargetsNormal, mLocusTargetNormal);
        makeAdapters(mMsiTargetsNormal, mLocusTargetNormal);
        evaluate(mOptions->normal, mMsiTargetsNormal, mLocusTargetNormal);
    }
    reportText();
    reportJSON();
    reportHTML();
}

void VisualMSI::evaluate(string filename, vector<MsiTarget*>& msiTargets, map<string, map<int, MsiTarget*>>& locusTarget){
    cerr << endl << "processing " << filename << endl;
    samFile *in;
    in = sam_open(filename.c_str(), "r");
    if (!in) {
        cerr << "ERROR: failed to open " << filename << endl;
        exit(-1);
    }

    bam_hdr_t* mBamHeader = sam_hdr_read(in);
    mOptions->bamHeader = mBamHeader;
    if (mBamHeader == NULL || mBamHeader->n_targets == 0) {
        cerr << "ERROR: this SAM file has no header " << filename << endl;
        exit(-1);
    }
    //BamUtil::dumpHeader(mBamHeader);

    bam1_t *b = NULL;
    b = bam_init1();
    int r;
    int count = 0;
    int lastTid = -1;
    int lastPos = -1;
    string chr;
    vector<int> locuses;
    while ((r = sam_read1(in, mBamHeader, b)) >= 0) {

        // unmapped reads, we just continue
        if(b->core.tid < 0 || b->core.pos < 0) {
            continue;
        }

        // check whether the BAM is sorted
        if(b->core.tid <lastTid || (b->core.tid == lastTid && b->core.pos <lastPos)) {
            // skip the -1:-1, which means unmapped
            if(b->core.tid >=0 && b->core.pos >= 0) {
                cerr << "ERROR: the input is unsorted. Found unsorted read in " << b->core.tid << ":" << b->core.pos << endl;
                cerr << "Please sort the input first." << endl << endl;
                exit(-1);
            }
        }

        if(b->core.tid != lastTid) {
            chr = string(mOptions->bamHeader->target_name[b->core.tid]);
            locuses.clear();
            if(locusTarget.count(chr)>0) {
                map<int, MsiTarget*>::iterator iter;
                for(iter = locusTarget[chr].begin(); iter != locusTarget[chr].end(); iter++)
                    locuses.push_back(iter->first);
            }
        }

        // for testing, we only process to some contig
        if(mOptions->maxContig>0 && b->core.tid>=mOptions->maxContig){
            b = bam_init1();
            break;
        }

        // if debug flag is enabled, show which contig we are start to process
        if(mOptions->debug && b->core.tid > lastTid) {
            cerr << "Starting contig " << b->core.tid << endl;
        }

        lastTid = b->core.tid;
        lastPos = b->core.pos;

        // for secondary alignments, we just skip it
        if(!BamUtil::isPrimary(b)) {
            continue;
        }

        // try to find whether this read is close to a MSI locus
        int msipos = -1;
        for(int l=0; l<locuses.size(); l++) {
            if(abs(b->core.isize) > 1000 || b->core.isize==0)
                continue;
            int bar = locuses[l];

            const int margin = 20;

            // read1
            if(b->core.pos < bar-margin && b->core.pos + b->core.isize > bar+margin) {
                msipos = bar;
                break;
            }

            // read2
            if(b->core.isize < 0 ) {
                if(b->core.mpos < bar-margin && b->core.mpos - b->core.isize > bar+margin) {
                    msipos = bar;
                    break;
                }
            }
        }

        if(msipos > 0) {
            addToMsiLocus(b, locusTarget[chr][msipos]);
            b = bam_init1();
        }
    }

    bam_destroy1(b);
    sam_close(in);

    stat(msiTargets, locusTarget);
}

void VisualMSI::addToMsiLocus(bam1_t* b, MsiTarget* t) {
    t->addRead(b);
}