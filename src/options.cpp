#include "options.h"
#include "util.h"

Options::Options(){
    input = "";
    targetFile = "";
    bamHeader = NULL;
    refFile = "";
    maxContig = 0;
    adapterLen = 12;
    targetInsertedSize = 100;
    depthReq = 10;
}

bool Options::validate() {
    if(input.empty()) {
        error_exit("input should be specified by --in1");
    } else {
        check_file_valid(input);
    }

    check_file_valid(refFile);
    check_file_valid(targetFile);

    if(!normal.empty())
        check_file_valid(normal);

    if(ends_with(refFile, ".gz") || ends_with(refFile, ".gz")) {
        cerr << "reference fasta file should not be compressed.\nplease unzip "<<refFile<<" and try again."<<endl;
        exit(-1);
    }

    if(adapterLen < 5 || adapterLen > 30)
        error_exit("Adapter length (--adapter_len) should be between 5 ~ 30.");

    if(depthReq < 1 || depthReq > 1000)
        error_exit("Depth requirement (--depth_req) should be between 1 ~ 1000.");

    if(targetInsertedSize < 20 || targetInsertedSize > 200)
        error_exit("Target inserted length (--target_inserted_len) should be between 20 ~ 200.");

    return true;
}