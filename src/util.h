#ifndef UTIL_H
#define UTIL_H

#include <stdlib.h>
#include <string>
#include <iostream>
#include <vector>
#include <sys/stat.h>
#include <algorithm>

using namespace std;

inline char complement(char base) {
    switch(base){
        case 'A':
        case 'a':
            return 'T';
        case 'T':
        case 't':
            return 'A';
        case 'C':
        case 'c':
            return 'G';
        case 'G':
        case 'g':
            return 'C';
        default:
            return 'N';
    }
}

inline bool starts_with( string const & value,  string const & starting)
{
    if (starting.size() > value.size()) return false;
    return  equal(starting.begin(), starting.end(), value.begin());
}

inline bool ends_with( string const & value,  string const & ending)
{
	if (ending.size() > value.size()) return false;
	return  equal(ending.rbegin(), ending.rend(), value.rbegin());
}

inline string trim(const string& str)
{
    string::size_type pos = str.find_first_not_of(' ');
    if (pos == string::npos)
    {
        return string("");
    }
    string::size_type pos2 = str.find_last_not_of(' ');
    if (pos2 != string::npos)
    {
        return str.substr(pos, pos2 - pos + 1);
    }
    return str.substr(pos);
}

inline int split(const string& str, vector<string>& ret_, string sep = ",")
{
    if (str.empty())
    {
        return 0;
    }

    string tmp;
    string::size_type pos_begin = str.find_first_not_of(sep);
    string::size_type comma_pos = 0;

    while (pos_begin != string::npos)
    {
        comma_pos = str.find(sep, pos_begin);
        if (comma_pos != string::npos)
        {
            tmp = str.substr(pos_begin, comma_pos - pos_begin);
            pos_begin = comma_pos + sep.length();
        }
        else
        {
            tmp = str.substr(pos_begin);
            pos_begin = comma_pos;
        }

        ret_.push_back(tmp);
        tmp.clear();
    }
    return 0;
}

inline string replace(const string& str, const string& src, const string& dest)
{
    string ret;

    string::size_type pos_begin = 0;
    string::size_type pos       = str.find(src);
    while (pos != string::npos)
    {
        ret.append(str.data() + pos_begin, pos - pos_begin);
        ret += dest;
        pos_begin = pos + 1;
        pos       = str.find(src, pos_begin);
    }
    if (pos_begin < str.length())
    {
        ret.append(str.begin() + pos_begin, str.end());
    }
    return ret;
}

inline string basename(const string& filename){
    string::size_type pos = filename.find_last_of('/');
    if (pos == string::npos)
        return filename;
    else if(pos == filename.length()-1)
        return ""; // a bad filename
    else
        return filename.substr(pos+1, filename.length() - pos - 1);
}

inline string dirname(const string& filename){
    string::size_type pos = filename.find_last_of('/');
    if (pos == string::npos) {
        return "./";
    } else
        return filename.substr(0, pos+1);
}

inline string joinpath(const string& dirname, const string& basename){
    if(dirname[dirname.length()-1] == '/'){
        return dirname + basename;
    } else {
        return dirname + "/" + basename;
    }
}

//Check if a string is a file or directory
inline bool file_exists(const  string& s)
{
    bool exists = false;
    if(s.length() > 0) {
        struct stat status;
        int result = stat( s.c_str(), &status );
        if(result == 0) {
            exists = true;
        }
    }
    return exists;
}


// check if a string is a directory
inline bool is_directory(const  string& path)
{
    bool isdir = false;
    struct stat status;
    // visual studion use _S_IFDIR instead of S_IFDIR
    // http://msdn.microsoft.com/en-us/library/14h5k7ff.aspx
#ifdef _MSC_VER
#define S_IFDIR _S_IFDIR
#endif
    stat( path.c_str(), &status );
    if ( status.st_mode &  S_IFDIR  ) {
        isdir = true;
    }
// #endif
    return isdir;
}

inline void check_file_valid(const  string& s) {
    if(!file_exists(s)){
        cerr << "ERROR: file '" << s << "' doesn't exist, quit now" << endl;
        exit(-1);
    }
    if(is_directory(s)){
        cerr << "ERROR: '" << s << "' is a folder, not a file, quit now" << endl;
        exit(-1);
    }
}

// Remove non alphabetic characters from a string
inline  string str_keep_alpha(const  string& s)
{
     string new_str;
    for( size_t it =0; it < s.size(); it++) {
        if(  isalpha(s[it]) ) {
            new_str += s[it];
        }
    }
    return new_str;
}


// Remove invalid sequence characters from a string
inline void str_keep_valid_sequence(  string& s, bool forceUpperCase = false)
{
    size_t total = 0;
    const char case_gap = 'a' - 'A';
    for( size_t it =0; it < s.size(); it++) {
        char c = s[it];
        if(forceUpperCase && c>='a' && c<='z') {
            c -= case_gap;
        }
        if(  isalpha(c) || c == '-' || c == '*' ) {
            s[total] = c;
            total ++;
        }
    }

    s.resize(total);
}

inline int find_with_right_pos(const string& str, const string& pattern, int start=0) {
    int pos = str.find(pattern, start);
    if (pos < 0)
        return -1;
    else
        return pos + pattern.length();
}

inline void str2upper(string& s){
    transform(s.begin(), s.end(), s.begin(), (int (*)(int))toupper);
}

inline void str2lower(string& s){
    transform(s.begin(), s.end(), s.begin(), (int (*)(int))tolower);
}

inline int hamming(const string& str1, const string& str2) {
    int diff = 0;
    int len1 = str1.length();
    int len2 = str2.length();
    for(int i=0; i<len1 && i<len2; i++) {
        if(str1[i] != str2[i])
            diff++;
    }
    diff += abs(len1 - len2);
    return diff;
}

inline char num2qual(int num) {
    if(num > 127 - 33)
        num = 127 - 33;
    if(num < 0)
        num = 0;

    char c = num + 33;
    return c;
}

inline void error_exit(const string& msg) {
    cerr << "ERROR: " << msg << endl;
    exit(-1);
}

/*
return the edit distance for string str1 and str2. Both str1 and str2 have a length of len
the gapPos will be returned to indicate whether there is an gap happens.
if gapPos > 0, there is a gap in str1, at the position gapPos
if gapPos < 0, there is a gap in str1, at the position gapPos
if gapPos = 0, there is no gap
the gap cannot happen at pos 0 or len-1
*/
inline int diffWithOneGap(const char* str1, const char* str2, unsigned int len, int& gapPos) {
    gapPos = 0;
    int diff = 0;
    for(unsigned int i=0; i<len; i++) {
        if(str1[i] != str2[i])
            diff++;
    }
    if(diff <= 1)
        return diff;

    int minDiff = 0x7FFFFFFF;
    int minDiffgapPos = 0;
    for(int gap = 1; gap < len-1; gap++) {
        int diff1 = 0;
        int diff2 = 0;
        for(unsigned int i=0; i<gap; i++) {
            if(str1[i] != str2[i]) {
                diff1++;
                diff2++;
            }
        }
        for(unsigned int i=gap; i<len-1; i++) {
            // gap at str1
            if(str1[i] != str2[i+1])
                diff1++;
            // gap at str2
            if(str1[i+1] != str2[i])
                diff2++;
        }
        if(diff1 < minDiff) {
            minDiff = diff1;
            minDiffgapPos = gap;
        }
        if(diff2 < minDiff) {
            minDiff = diff2;
            minDiffgapPos = -gap;
        }
    }

    if(minDiff < diff - 1) {
        gapPos = minDiffgapPos;
        return minDiff;
    }

    return diff;
}

inline bool diffWithOneGapTest() {
    string str1 = "ATCGATCGATCGATCG";
    string str2 = "ATCGATCATCGATCGA";
    const char* data1 = str1.c_str();
    const char* data2 = str2.c_str();
    int diff, gapPos;
    diff = diffWithOneGap(data1, data2, 16, gapPos);
    cerr << diff << ":" << gapPos << endl;
    if(diff != 0 || gapPos != -7)
        return false;

    diff = diffWithOneGap(data1, data2, 6, gapPos);
    cerr << diff << ":" << gapPos << endl;
    if(diff != 0 || gapPos != 0)
        return false;

    diff = diffWithOneGap(data1, data2, 8, gapPos);
    cerr << diff << ":" << gapPos << endl;
    if(diff != 1 || gapPos != 0)
        return false;

    str1 = "ATCGACGTCGAA";
    str2 = "ATCGTACGTCGA";

    diff = diffWithOneGap(data1, data2, 12, gapPos);
    cerr << diff << ":" << gapPos << endl;
    if(diff != 0 || gapPos != 4)
        return false;

    diff = diffWithOneGap(data1+3, data2+4, 5, gapPos);
    cerr << diff << ":" << gapPos << endl;
    if(diff != 1 || gapPos != 0)
        return false;

    return true;
}

#endif /* UTIL_H */
