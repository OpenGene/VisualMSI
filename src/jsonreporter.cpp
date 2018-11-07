#include "jsonreporter.h"

JsonReporter::JsonReporter(Options* opt){
    mOptions = opt;
}

JsonReporter::~JsonReporter(){
}

void JsonReporter::reportMsiTarget(ofstream& ofs, MsiTarget* tumor, MsiTarget* normal) {
    ofs << "\t\t\"" << tumor->mName << "\":{" << endl;

    if(normal)
        ofs << "\t\t\t\"emd_distance\":" << tumor->emdWith(normal) << "," << endl;

    ofs << "\t\t\t\"tumor_entropy\":" << tumor->entropy() << "," << endl;
    if(normal)
        ofs << "\t\t\t\"normal_entropy\":" << normal->entropy() << "," << endl;

    ofs << "\t\t\t\"tumor_supporting_reads\":" << tumor->mSupportingReads << "," << endl;
    if(normal)
        ofs << "\t\t\t\"normal_supporting_reads\":" << normal->mSupportingReads << "," << endl; 

    if(tumor->mSupportingReads < mOptions->depthReq || (normal && normal->mSupportingReads < mOptions->depthReq))
        ofs << "\t\t\t\"qc_result\":\"failed\"," << endl; 
    else
        ofs << "\t\t\t\"qc_result\":\"passed\"," << endl; 

    ofs << "\t\t\t\"tumor_tlen_histogram\":[";
    for(int i=1; i<MAX_HIST_BIN; i++) {
        ofs << tumor->mHistogram[i];
        if(i!=MAX_HIST_BIN - 1)
            ofs << ",";
    }
    ofs <<"]," << endl;

    if(normal) {
        ofs << "\t\t\t\"normal_tlen_histogram\":[";
        for(int i=1; i<MAX_HIST_BIN; i++) {
            ofs << normal->mHistogram[i];
            if(i!=MAX_HIST_BIN - 1)
                ofs << ",";
        }
        ofs <<"]," << endl;
    }

    ofs << "\t\t\t\"start\":" << tumor->mStart << "," << endl;
    ofs << "\t\t\t\"end\":" << tumor->mEnd << "," << endl;
    ofs << "\t\t\t\"left_adapter\":\"" << tumor->mLeftAdapter << "\"," << endl;
    ofs << "\t\t\t\"right_adapter\":\"" << tumor->mRightAdapter << "\"," << endl;
    ofs << "\t\t\t\"inserted_on_ref\":\"" << tumor->mInserted << "\"" << endl;
    ofs << "\t\t}";
}

extern string command;
void JsonReporter::report(vector<MsiTarget*>& msiTargetsTumor, vector<MsiTarget*>& msiTargetsNormal){
    bool hasNormal = !mOptions->normal.empty() && msiTargetsNormal.size() == msiTargetsTumor.size();

	ofstream ofs;
    ofs.open(mOptions->jsonFile, ifstream::out);
    ofs << "{" << endl;

    ofs << "\t" << "\"summary\": {" << endl;
    ofs << "\t\t\"tumor_file\":\"" << mOptions->input << "\"," << endl;
    if(hasNormal)
        ofs << "\t\t\"normal_file\":\"" << mOptions->normal << "\"," << endl;
    ofs << "\t\t\"target_file\":\"" << mOptions->targetFile << "\"," << endl;
    ofs << "\t\t\"total_sites\":" << msiTargetsTumor.size();
    ofs << endl;
    ofs << "\t" << "}," << endl;

    ofs << "\t" << "\"detail\": {" << endl;

    for(int i=0; i<msiTargetsTumor.size();i++){
        if(hasNormal)
            reportMsiTarget(ofs, msiTargetsTumor[i], msiTargetsNormal[i]);
        else
            reportMsiTarget(ofs, msiTargetsTumor[i], NULL);
    	if(i != msiTargetsTumor.size() - 1)
    		ofs << "," << endl;
    }

    ofs << endl << "\t" << "}," << endl;

    ofs << "\t\"command\": " << "\"" << command << "\"" << endl;

    ofs << "}";

    ofs.close();
}