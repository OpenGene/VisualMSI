#ifndef JSON_REPORTER_H
#define JSON_REPORTER_H

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include "options.h"
#include <fstream>
#include "msitarget.h"

using namespace std;

class JsonReporter{
public:
    JsonReporter(Options* opt);
    ~JsonReporter();

    void report(vector<MsiTarget*>& msiTargetsTumor, vector<MsiTarget*>& msiTargetsNormal);
    void reportMsiTarget(ofstream& ofs, MsiTarget* tumor, MsiTarget* normal);

private:
    Options* mOptions;
};


#endif