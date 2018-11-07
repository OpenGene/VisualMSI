#ifndef PAIR_H
#define PAIR_H

#include <stdio.h>
#include <stdlib.h>
#include "util.h"
#include "htslib/sam.h"
#include "options.h"

using namespace std;

class Pair {
public:
    Pair();
    ~Pair();

    enum MapType{Unknown, ProperlyMapped, CrossRefMapped, OnlyLeftMapped, OnlyRightMapped, NoneMapped};

    int getLeftRef();
    int getRightRef();
    int getLeftPos();
    int getRightPos();
    int getTLEN();

    void setLeft(bam1_t *b);
    void setRight(bam1_t *b);
    bool pairFound();
    MapType getMapType();
    string getUMI();
    string getQName();
    string getLeftCigar();
    string getRightCigar();

    string getInsertedByPCR(string leftAdapter, string rightAdapter);
    static string getInsertedByPCR(string seq, string leftAdapter, string rightAdapter);
    
    void dump();

    static bool test();


public:
    bam1_t *mLeft;
    bam1_t *mRight;

private:
    int mTLEN;
    MapType mMapType;
    string mUMI;
    string mLeftCigar;
    string mRightCigar;
    Options* mOptions;
};

#endif