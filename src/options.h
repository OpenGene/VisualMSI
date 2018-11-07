#ifndef OPTIONS_H
#define OPTIONS_H

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>
#include <map>
#include "htslib/sam.h"

using namespace std;


class Options{
public:
    Options();
    bool validate();

public:
    bam_hdr_t* bamHeader;
    string input;
    string normal;
    string targetFile;
    string refFile;
    int maxContig;
    bool debug;
    // json file
    string jsonFile;
    // html file
    string htmlFile;

    // thresholds
    int adapterLen;
    int targetInsertedSize;
    int depthReq;
};

#endif