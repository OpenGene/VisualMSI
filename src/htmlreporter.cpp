#include "htmlreporter.h"
#include <sstream>
#include <chrono>
#include "common.h"

HtmlReporter::HtmlReporter(Options* opt){
    mOptions = opt;
}

HtmlReporter::~HtmlReporter(){
}

extern string command;

void HtmlReporter::printHeader(ofstream& ofs){
    ofs << "<html><head><meta http-equiv=\"content-type\" content=\"text/html;charset=utf-8\" />";
    ofs << "<title>VisualMSI report, " + getCurrentSystemTime() + " </title>";
    printJS(ofs);
    printCSS(ofs);
    ofs << "</head>";
    ofs << "<body><div id='container'>";
    ofs << "<div class='header'><a href='https://github.com/OpenGene/VisualMSI' target='_blank'>VisualMSI</a> report</div>\n";
}

void HtmlReporter::printCSS(ofstream& ofs){
    ofs << "<style type=\"text/css\">" << endl;
    ofs << "td {border:1px solid #dddddd;padding:5px;font-size:12px;}" << endl;
    ofs << "table {border:1px solid #999999;padding:2x;border-collapse:collapse; width:800px}" << endl;
    ofs << ".col1 {width:240px; font-weight:bold;}" << endl;
    ofs << ".adapter_col {width:500px; font-size:10px;}" << endl;
    ofs << "img {padding:30px;}" << endl;
    ofs << "#menu {font-family:Consolas, 'Liberation Mono', Menlo, Courier, monospace;}" << endl;
    ofs << "#menu a {color:#0366d6; font-size:18px;font-weight:600;line-height:28px;text-decoration:none;font-family:-apple-system, BlinkMacSystemFont, 'Segoe UI', Helvetica, Arial, sans-serif, 'Apple Color Emoji', 'Segoe UI Emoji', 'Segoe UI Symbol'}" << endl;
    ofs << "a:visited {color: #999999}" << endl;
    ofs << ".alignleft {text-align:left;}" << endl;
    ofs << ".alignright {text-align:right;}" << endl;
    ofs << ".figure {width:800px;height:500px;}" << endl;
    ofs << ".header {color:#ffffff;padding:5px;height:20px;background:#663355;}" << endl;
    ofs << ".section_title {color:#ffffff;font-size:18px;padding:5px;text-align:left;background:#663355; margin-top:10px;}" << endl;
    ofs << ".subsection_title {font-size:16px;padding:5px;margin-top:10px;text-align:left;color:#663355}" << endl;
    ofs << "#container {text-align:left;padding:3px 3px 3px 10px;font-family:Arail,'Liberation Mono', Menlo, Courier, monospace;}" << endl;
    ofs << ".menu_item {text-align:left;padding-top:5px;font-size:18px;}" << endl;
    ofs << ".highlight {text-align:left;padding-top:30px;padding-bottom:30px;font-size:20px;line-height:35px;}" << endl;
    ofs << "#helper {text-align:left;border:1px dotted #fafafa;color:#777777;font-size:12px;}" << endl;
    ofs << "#footer {text-align:left;padding:15px;color:#ffffff;font-size:10px;background:#663355;font-family:Arail,'Liberation Mono', Menlo, Courier, monospace;}" << endl;
    ofs << ".kmer_table {text-align:center;font-size:8px;padding:2px;}" << endl;
    ofs << ".kmer_table td{text-align:center;font-size:8px;padding:0px;color:#ffffff}" << endl;
    ofs << ".sub_section_tips {color:#666666;font-size:8px;padding-left:5px;padding-bottom:3px;}" << endl;
    ofs << "</style>" << endl;
}

void HtmlReporter::printJS(ofstream& ofs){
    ofs << "<script src='https://cdn.plot.ly/plotly-latest.min.js'></script>" << endl;
    ofs << "\n<script type=\"text/javascript\">" << endl;
    ofs << "    function showOrHide(divname) {" << endl;
    ofs << "        div = document.getElementById(divname);" << endl;
    ofs << "        if(div.style.display == 'none')" << endl;
    ofs << "            div.style.display = 'block';" << endl;
    ofs << "        else" << endl;
    ofs << "            div.style.display = 'none';" << endl;
    ofs << "    }" << endl;
    ofs << "</script>" << endl;
}

const string HtmlReporter::getCurrentSystemTime()
{
  auto tt = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
  struct tm* ptm = localtime(&tt);
  char date[60] = {0};
  sprintf(date, "%d-%02d-%02d      %02d:%02d:%02d",
    (int)ptm->tm_year + 1900,(int)ptm->tm_mon + 1,(int)ptm->tm_mday,
    (int)ptm->tm_hour,(int)ptm->tm_min,(int)ptm->tm_sec);
  return std::string(date);
}

void HtmlReporter::printFooter(ofstream& ofs){
    ofs << "\n</div>" << endl;
    ofs << "<div id='footer'> ";
    ofs << "<p>"<<command<<"</p>";
    ofs << "VisualMSI " << VERSION_NUMBER << ", at " << getCurrentSystemTime() << " </div>";
    ofs << "</body></html>";
}

template <class T>
string HtmlReporter::list2string(T* list, int size) {
    stringstream ss;
    for(int i=0; i<size; i++) {
        ss << list[i];
        if(i < size-1)
            ss << ",";
    }
    return ss.str();
}

void HtmlReporter::reportMsiTarget(ofstream& ofs, MsiTarget* tumor, MsiTarget* normal) {
    if(tumor->mSupportingReads < mOptions->depthReq) {
        ofs << "<div class='subsection_title'> " << tumor->mName << ": tumor data depth is too low, require " << mOptions->depthReq;
        ofs << ", but we have only " << tumor->mSupportingReads << "</div>\n";
        return;
    } else if(normal && normal->mSupportingReads < mOptions->depthReq) {
        ofs << "<div class='subsection_title'> " << normal->mName << ": normal data depth is too low, require " << mOptions->depthReq;
        ofs << ", but have only " << normal->mSupportingReads << "</div>\n";
        return;
    }
    // quality
    string subsection = tumor->mName + "(";
    if(normal)
        subsection += "tumor entropy: " + to_string(tumor->entropy()) + ", normal entropy: " + to_string(normal->entropy()) + ", EMD distance: " + to_string(tumor->emdWith(normal)) + ")";
    else
        subsection += " entropy: " + to_string(tumor->entropy()) + ")"; 

    string divName = replace(subsection, " ", "_");
    divName = replace(divName, ":", "_");
    string title = tumor->mName;
    string seq = tumor->mLeftAdapter + " " + tumor->mInserted + " " + tumor->mRightAdapter;

    ofs << "<div class='subsection_title'><a title='click to hide/show' onclick=showOrHide('" << divName << "')>" + subsection + "</a></div>\n";
    ofs << "<div id='" + divName + "'>\n";
    ofs << "<div class='sub_section_tips'>" << tumor->mChr << ":" << tumor->mStart << "-" << tumor->mEnd << "</div>\n";
    ofs << "<div class='sub_section_tips'>Left adapter: " << tumor->mLeftAdapter << "</div>\n";
    ofs << "<div class='sub_section_tips'>Right adapter: " << tumor->mRightAdapter << "</div>\n";
    ofs << "<div class='sub_section_tips'>Reference: ..." << seq << "...</div>\n";
    ofs << "<div class='figure' id='plot_" + divName + "'></div>\n";
    ofs << "</div>\n";

    int minBin = tumor->mMinBin;
    if(normal)
        minBin = min(minBin, tumor->mMinBin);

    int maxBin = tumor->mMaxBin;
    if(normal)
        maxBin = max(maxBin, tumor->mMaxBin);
    int bins = maxBin-minBin+1;

    if(bins<20) {
        minBin -= (20-bins)/2;
        bins = 20;
    }

    int *x = new int[bins];
    for(int i=0; i<bins; i++)
        x[i] = minBin + i;
    
    string alphabets[2] = {"tumor", "normal"};
    string colors[2] = {"rgba(255,0, 0,1.0)", "rgba(0, 0, 255,1.0)"};
    ofs << "\n<script type=\"text/javascript\">" << endl;
    string json_str = "var data=[";

    int channel = 1;

    if(normal) {
        channel = 2;
    }

    double tumorTotal = 0.0;
    double normalTotal = 0.0;
    double* tumorHistogram = new double[bins];
    double* normalHistogram = new double[bins];

    // normalization of the graph
    for(int i=0; i<bins; i++)
        tumorTotal += tumor->mHistogram[i+minBin];
    for(int i=0; i<bins; i++)
        tumorHistogram[i] = tumor->mHistogram[i+minBin] / tumorTotal;

    if(normal) {
        for(int i=0; i<bins; i++)
            normalTotal += normal->mHistogram[i+minBin];
        for(int i=0; i<bins; i++)
            normalHistogram[i] = normal->mHistogram[i+minBin] / normalTotal;
    }

    // four bases
    for (int b = 0; b<channel; b++) {
        string base = alphabets[b];
        json_str += "{";
        json_str += "x:[" + list2string(x, bins) + "],";
        if(b == 0)
            json_str += "y:[" + list2string(tumorHistogram, bins) + "],";
        else
            json_str += "y:[" + list2string(normalHistogram, bins) + "],";
        json_str += "name: '" + base + "',";
        json_str += "fill: 'tozeroy',";
        json_str += "line:{color:'" + colors[b] + "', width:0}\n";
        json_str += "},";
    }
    json_str += "];\n";
    json_str += "var layout={title:'" + title + "', xaxis:{title:'inserted size by PCR simulation'";
    json_str += "}, yaxis:{title:'ratio'}};\n";
    json_str += "Plotly.newPlot('plot_" + divName + "', data, layout);\n";

    ofs << json_str;
    ofs << "</script>" << endl;

    delete[] x;
}

void HtmlReporter::report(vector<MsiTarget*>& msiTargetsTumor, vector<MsiTarget*>& msiTargetsNormal){
	ofstream ofs;
    ofs.open(mOptions->htmlFile, ifstream::out);
    printHeader(ofs);

    bool hasNormal = !mOptions->normal.empty() && msiTargetsNormal.size() == msiTargetsTumor.size();

    for(int i=0; i<msiTargetsTumor.size();i++){
        if(hasNormal)
            reportMsiTarget(ofs, msiTargetsTumor[i], msiTargetsNormal[i]);
        else
            reportMsiTarget(ofs, msiTargetsTumor[i], NULL);
    }

    printFooter(ofs);

    ofs.close();
}