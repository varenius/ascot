#ifndef AUXFUNC_H
#define AUXFUNC_H

#include<iostream>
#include<vector>
#include<iterator>
#include<algorithm>
#include<string>
#include<fstream>
#include<sstream>
#include<stdexcept>
#include<map>
#include "matrix.h"
#include<boost/algorithm/string.hpp>
#include<boost/math/special_functions/modf.hpp>
#include<boost/program_options.hpp>

#include <cstdlib>
#include <libconfig.h++>

using namespace libconfig;

using namespace std;

bool gt( double i, double j );
bool lt( double i, double j );
bool ge( double i, double j );
bool le( double i, double j );
bool eq( double i, double j );
bool ne( double i, double j );

char* cStr( string str );

bool is_file_exist( const char *fileName );

void replace_string_in_place( std::string& subject, const std::string& search,
                              const std::string& replace );

string D2E(string str);

string remove_spaces(string str);

string remove_spaces_end(const string s);

Setting& get_list_element(Setting &setting, string key, int argument=0);

std::vector<std::pair<std::string, std::string> > get_baselines(Setting &setting);

bool includes_baseline(std::vector<std::pair<std::string, std::string> > bl, std::string sta1, std::string sta2 );

bool file_exists (const std::string& name);

bool directory_exists (const std::string& name);

/**
*  \b Description: \n
*        Function saving local filelist in a vector
*  \param [in] [std::string] path e.g. "/home/ascot/"
*  \return [vector<std::string>] all entries in remote location
*/
std::vector<std::string> list_local_dir(std::string dir);

// set read/write for user and group and read-only for others
void chmod_urw_grw_or( std::string file );

// set executable for all (user+group+other)
void chmod_ax( std::string file );

bool ends_with(const string& str, const string& ending );

vector<string> get_tokens(const string &line);

// get stations from groups, selected in config file
void group2names( Setting &params, Setting &groups,
                  std::string selected_group_name,
                  std::vector< std::string > &names, 
                  std::vector< std::string > &param_lst );
// string to double, works in every case
// solves the problem using QApllication and std::stod
double s2d(string str);


double azimuth_diff(double az1, double az2);

// ---------------------------------------------------------------------------
// functions with template classes have to be implemented in the header file


// split a string in its elements (sperated by spaces)
// input string s should contain only data of one type, i.e., string, int,
// double ... otherwise the result in out might be corrupted
void splitStringSpaces2( const string &s, vector<int> &out );
template <class T>
void splitStringSpaces( const string &s, vector<T> &out )
{
    istringstream iss(s);
    T sub;

    out.erase( out.begin(), out.end() );

    while( iss >> sub )
        out.push_back( sub );
}

// ---------------------------------------------------------------------------
template <typename T>
void show_vector(vector<T> vec) {
    copy ( vec.begin () , vec.end () , ostream_iterator <T>( cerr , "|" ) ); cerr << endl;
}  
// ---------------------------------------------------------------------------
template <typename T>
void remove_duplicates(std::vector<T>& vec)
{
  std::sort(vec.begin(), vec.end());
  vec.erase(std::unique(vec.begin(), vec.end()), vec.end());
}
// ---------------------------------------------------------------------------
#endif  // AUXFUNC_H
