#include "reference.h"
#include "util.h"

Reference* Reference::mInstance = NULL;

Reference* Reference::instance(Options* opt) {
    if(mInstance == NULL)
        mInstance = new Reference(opt);
    
    return mInstance;
}

Reference::Reference(Options* opt) {
    mOptions = opt;
    mRef = NULL;
    if(!mOptions->refFile.empty()) {
        mRef = new FastaReader(mOptions, mOptions->refFile);
        mRef->readAll();
    }
    mLastBamContig = -1;
    mLastData = NULL;
    mLastLen = -1;
}

Reference::~Reference() {
    if(mRef) {
        delete mRef;
        mRef = NULL;
    }
    mInstance = NULL;
}

string Reference::getContig(string contigName) {
    if(mRef->mAllContigs.count(contigName) == 0)
        return string();
    else
        return mRef->mAllContigs[contigName];
}

const char* Reference::getData(int bamContig, int pos, int len) {
    if(mRef == NULL)
        return NULL;
    if(mOptions->bamHeader == NULL)
        return NULL;

    if(mLastBamContig == bamContig && mLastData!=NULL) {
        if(pos + len >= mLastLen)
            return NULL;
        else
            return mLastData + pos;
    }

    // get contig name from bam header
    string contigName(mOptions->bamHeader->target_name[bamContig]);

    mLastBamContig = bamContig;

    if(mRef->mAllContigs.count(contigName) == 0) {
        static bool reported = false;
        if(!reported)
            cerr << "contig " << contigName << " not found in the reference, please make sure your reference is correct" << endl;
        reported = true;
        mLastData = NULL;
        return NULL;
    }

    if(pos + len >= mRef->mAllContigs[contigName].length()){
        static bool reported = false;
        if(!reported)
            cerr << "contig " << contigName << " doesn't match the length in the reference, please make sure your reference is correct" << endl;
        mLastData = NULL;
        return NULL;
    }

    mLastData = mRef->mAllContigs[contigName].c_str();
    mLastLen = mRef->mAllContigs[contigName].length();
    return mLastData + pos;
}