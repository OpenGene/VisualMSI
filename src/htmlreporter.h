#ifndef HTML_REPORTER_H
#define HTML_REPORTER_H

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include "options.h"
#include <fstream>
#include "msitarget.h"

using namespace std;

class HtmlReporter{
public:
    HtmlReporter(Options* opt);
    ~HtmlReporter();

    void report(vector<MsiTarget*>& msiTargetsTumor, vector<MsiTarget*>& msiTargetsNormal);
    void reportMsiTarget(ofstream& ofs, MsiTarget* tumor, MsiTarget* normal);

private:
    const string getCurrentSystemTime();
    void printHeader(ofstream& ofs);
    void printCSS(ofstream& ofs);
    void printJS(ofstream& ofs);
    void printFooter(ofstream& ofs);
    template <class T>
    string list2string(T* list, int size);

private:
    Options* mOptions;
};


#endif