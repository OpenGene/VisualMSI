#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "cmdline.h"
#include "common.h"
#include <sstream>
#include "util.h"
#include "visualmsi.h"
#include "options.h"
#include "reference.h"
#include "unittest.h"

using namespace std;

string command;

int main(int argc, char* argv[]){
    if (argc == 2 && strcmp(argv[1], "test")==0){
        UnitTest tester;
        tester.run();
        return 0;
    }

    if (argc == 2 && (strcmp(argv[1], "-v")==0 || strcmp(argv[1], "--version")==0)){
        cerr << "visualmsi " << VERSION_NUMBER << endl;
        return 0;
    }

    cmdline::parser cmd;
    // input/output
    cmd.add<string>("in", 'i', "input sorted bam/sam file for the case (tumor) sample. STDIN will be read from if it's not specified", false, "-");
    cmd.add<string>("normal", 'n', "input sorted bam/sam file for the paired normal sample (tumor-normal mode). If not specified, VisualMSI will run in case-only mode.", false, "");
    cmd.add<string>("target", 't', "the TSV file (chrom, position, name) to give the MSI targets", true, "");
    cmd.add<string>("ref", 'r', "reference fasta file name (should be an uncompressed .fa/.fasta file)", true, "");
    
    // UMI
    //cmd.add<string>("umi_prefix", 'u', "the prefix for UMI, if it has. None by default. Check the README for the defails of UMI formats.", false, "");

    cmd.add<int>("adapter_len", 'a', "set the length of the adapter for PCR simulation (5~30). Default 12 means the left and right adapter both have 12 bp.", false, 12);
    cmd.add<int>("target_inserted_len", 'l', "set the distance on reference of the two adapters for PCR simulation (20~200). Default 60 means: <left adapter><60 bp inserted><right adapter>", false, 60);
    cmd.add<int>("depth_req", 'd', "set the minimum depth requirement for each MSI locus (1~1000). Default 10 means 10 supporting reads/pairs are required.", false, 10);

    // reporting
    cmd.add<string>("json", 'j', "the json format report file name", false, "msi.json");
    cmd.add<string>("html", 'h', "the html format report file name", false, "msi.html");

    // debugging
    cmd.add("debug", 0, "output some debug information to STDERR.");

    cmd.parse_check(argc, argv);

    Options opt;
    opt.input = cmd.get<string>("in");
    opt.normal = cmd.get<string>("normal");
    opt.targetFile = cmd.get<string>("target");
    opt.refFile = cmd.get<string>("ref");
    opt.debug = cmd.exist("debug");
    opt.adapterLen = cmd.get<int>("adapter_len");
    opt.depthReq = cmd.get<int>("depth_req");
    opt.targetInsertedSize = cmd.get<int>("target_inserted_len");

    // reporting
    opt.jsonFile = cmd.get<string>("json");
    opt.htmlFile = cmd.get<string>("html");

    opt.validate();
    
    time_t t1 = time(NULL);

    // loading reference
    Reference* reference = NULL;
    if(!opt.refFile.empty()) {
        cerr << "loading reference data:" << endl;
        reference = Reference::instance(&opt);
    }

    stringstream ss;
    for(int i=0;i<argc;i++){
        ss << argv[i] << " ";
    }
    command = ss.str();

    VisualMSI vmsi(&opt);
    vmsi.run();

    if(reference) {
        delete reference;
        reference=NULL;
    }

    time_t t2 = time(NULL);
    cerr << endl << command << endl;
    cerr << "VisualMSI v" << VERSION_NUMBER << ", time used: " << (t2)-t1 << " seconds" << endl;

    return 0;
}