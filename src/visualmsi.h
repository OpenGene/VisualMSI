#ifndef VISUALMSI_H
#define VISUALMSI_H

#include <stdio.h>
#include <stdlib.h>
#include "util.h"
#include "htslib/sam.h"
#include "options.h"
#include "pair.h"
#include <map>
#include <vector>
#include "msitarget.h"

using namespace std;

class VisualMSI {
public:
    VisualMSI(Options *opt);
    ~VisualMSI();

    void run();
    void evaluate(string filename, vector<MsiTarget*>& msiTargets, map<string, map<int, MsiTarget*>>& locusTarget);

private:
    void addToMsiLocus(bam1_t* b, MsiTarget* t);
    void calcLocusTarget(vector<MsiTarget*>& msiTargets, map<string, map<int, MsiTarget*>>& locusTarget);
    void makeAdapters(vector<MsiTarget*>& msiTargets, map<string, map<int, MsiTarget*>>& locusTarget);
    void stat(vector<MsiTarget*>& msiTargets, map<string, map<int, MsiTarget*>>& locusTarget);
    void reportJSON();
    void reportHTML();
    void reportText();

private:
    Options *mOptions;
    vector<MsiTarget*> mMsiTargetsTumor;
    map<string, map<int, MsiTarget*>> mLocusTargetTumor;
    vector<MsiTarget*> mMsiTargetsNormal;
    map<string, map<int, MsiTarget*>> mLocusTargetNormal;
};

#endif