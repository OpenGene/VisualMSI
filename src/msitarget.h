#ifndef MSITARGET_H
#define MSITARGET_H

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <map>
#include "sequence.h"
#include "pair.h"
#include "htslib/sam.h"

#define MAX_HIST_BIN 256

using namespace std;

class MsiTarget{
public:
    MsiTarget(string chr, int start, int end, string name);
    ~MsiTarget();

    void addLen(int len);
    void addRead(bam1_t *b);
    static vector<MsiTarget*> parseBed(string filename);
    void print();
    void stat();
    double entropy();
    double emdWith(MsiTarget* other);

public:
    string mChr;
    int mStart;
    int mEnd;
    int mCenter;
    string mName;
    int* mHistogram;
    int mMaxBin;
    int mMinBin;
    int mSupportingReads;
    string mLeftAdapter;
    string mRightAdapter;
    string mInserted;
    map<string, Pair*> mPairs;
};


#endif
