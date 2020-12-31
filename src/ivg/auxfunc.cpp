#include "auxfunc.h"
#include <dirent.h>
#include "logger.h"

bool gt( double i, double j )
{
    return i > j;
}

bool lt( double i, double j )
{
    return i < j;
}

bool ge( double i, double j )
{
    return i >= j;
}

bool le( double i, double j )
{
    return i <= j;
}

bool eq( double i, double j )
{
    return i == j;
}

bool ne( double i, double j )
{
    return i != j;
}


char* cStr( string str )
{
    vector<char> buffer( str.begin(), str.end());
    char *key = &buffer[0];

    return key;
}

void replace_string_in_place( std::string& subject, const std::string& search,
                              const std::string& replace )
{
    size_t pos = 0;
    while ((pos = subject.find(search, pos)) != std::string::npos)
    {
        subject.replace(pos, search.length(), replace);
        pos += replace.length();
    }
}

bool is_file_exist(const char *fileName)
{
    std::ifstream infile(fileName);
    return infile.good();
}

string D2E(string str)
{
    size_t posi=str.find("D");
    if(posi!=string::npos)
    {
        str.replace(posi,1,"E");
    }

    return(str);
}


string remove_spaces(string str)
{
    string::size_type pos = 0;
    bool spacesLeft = true;

    while(spacesLeft)
    {
        pos = str.find(" ");
        if( pos != string::npos )
        {
            str.erase( pos, 1 );
        }
        else
        {
            spacesLeft = false;
        }
    }

    return(str);
}

string remove_spaces_end(const string s)
{
    int last = s.size() - 1;
    while (last >= 0 && s[last] == ' ')
        --last;
    return s.substr(0, last + 1);
}


Setting& get_list_element(Setting &setting, string key, int argument)
{

    //Detecting which type data for e.g. OCEAN LOADING
    int index;
    for(index=0; index<setting.getLength(); index++)
    {
        if( (const char *)setting[index][argument] == key)
        {
            Setting &defs = setting[index];
            return(defs);
        }
    }

}


std::vector<std::pair<std::string, std::string> > get_baselines(Setting &list)
{

    std::vector<std::pair<std::string, std::string>> bl (0);

    bl.resize(list.getLength());

    for(int i=0; i < list.getLength(); i++){  
        std::string first = list[i][0];
        std::string second = list[i][1];
        bl[i] = make_pair( first, second  );
    }
    return bl;
}

bool includes_baseline(std::vector<std::pair<std::string, std::string> > baselines, std::string sta1, std::string sta2){
   
    std::pair<std::string, std::string> ref1 = make_pair( sta1, sta2 );
    std::pair<std::string, std::string> ref2 = make_pair( sta2, sta1 );
    for( std::pair<std::string, std::string>& bl: baselines ){
        if( ref1 == bl || ref2 == bl){
            return true;
        } 
    }
    
    return false;
}

// ...........................................................................
double s2d(string str)
// ...........................................................................
{
    double d;
    stringstream stream(str);
    stream >> d;
    return d;
}

// ...........................................................................
double azimuth_diff(double az1, double az2){
// ...........................................................................
    double slew = 0.0;
    double alpha1 = az1 - az2;

    double alpha2 = -boost::math::sign(alpha1)*(360-abs(alpha1));

    // it is assumed the shorter angle is the correct one. Only not satisfied in pathological cases
    if( abs(alpha1) <= abs(alpha2) ){
        slew = alpha1;
    } else {
        slew = alpha2;
    }
    
    return slew;
}
// ...........................................................................
bool file_exists (const std::string& name) 
// ...........................................................................
{
  struct stat buffer;   
  return (stat (name.c_str(), &buffer) == 0); 
}

// ...........................................................................
bool directory_exists (const std::string& dir)
// ...........................................................................
{
    struct stat sb;
    return (stat(dir.c_str(), &sb) == 0 && S_ISDIR(sb.st_mode));
}

// .............................................................................
std::vector<std::string> list_local_dir(std::string dir){
// .............................................................................
    std::vector<std::string> ls;
    struct stat sb;
    
    //ckeck wether dir is an existing dir
    if (stat(dir.c_str(), &sb) == 0 && S_ISDIR(sb.st_mode)){
    
        struct dirent **namelist;
        int i,n;

        const char * c = dir.c_str();

        n = scandir(c, &namelist, 0, alphasort); //get all entries in the directory

        // check wether there are any entries
        if (n < 0)
            log<WARNING> ("!!! Scandir:") % dir % " is an empty directory";
        else 
        {   
            for (i = 0; i < n; i++)
            {
                std::string s = namelist[i]->d_name;
                ls.push_back(s);
                //std::cout << s << std::endl;
                free(namelist[i]);
            }
        }

       free(namelist);
    }
    else
        log<WARNING> ("!!! ") % dir % " is not a valid path";

   return ls;
}

// ...........................................................................
void chmod_urw_grw_or( std::string file )
// ...........................................................................
{
  int amode = S_IWUSR | S_IRUSR | S_IWGRP | S_IRGRP | S_IROTH;
  if ( chmod( file.c_str(),amode ) != 0 )
     cerr << "WARNING: permissions of " << file << " could not be changed" << endl;
}

// ...........................................................................
void chmod_ax( std::string file )
// ...........................................................................
{
  int amode = S_IXUSR | S_IXGRP | S_IXOTH;
  if ( chmod( file.c_str(),amode ) != 0 )
     cerr << "WARNING: permissions of " << file << " could not be changed" << endl;
}

// ...........................................................................
bool ends_with(const string& str,const string& ending )
// ...........................................................................
{
    if (str.length() >= ending.length())
    {
        std::string str_end = str.substr( str.length()-ending.length(), str.length() );
        return  (ending.compare( str_end ) == 0);
    }
    else
        return false;
}

// ...........................................................................
vector<string> get_tokens(const string &line)
// ...........................................................................
{
    string val;
    vector<string> tokens;
    stringstream tokenizer(line);
    while(tokenizer >> val)
        tokens.push_back(val);
    
    return tokens;
}
// ...........................................................................
void group2names( Setting &params, Setting &groups, 
                  std::string selected_group_name,
                  std::vector< std::string > &names, 
                  std::vector< std::string > &param_lst )
// ...........................................................................
{
   // loop over list of this config-file entry

   // group index 0  -> all stations, source or what ever
   //             >0 -> this group
   //             <0 -> all but this group
   //int name_idx = params[ selected_group_name ];
   int name_idx = params.lookup( selected_group_name );
   
   param_lst.clear();

   // select all param names of this class (stations, sources, ... ) within this session
   if( name_idx == 0)
      param_lst = names;
            
   // select only those param-names within the requested group
   else if( name_idx > 0 )
   {
       //for( int i=0;i<groups[selected_group_name][ name_idx-1 ].getLength(); ++i )
       for( int i=0;i<groups.lookup(selected_group_name)[ name_idx-1 ].getLength(); ++i )
       {
          //std::string elem = groups[selected_group_name][ name_idx-1 ][ i ];
          std::string elem = groups.lookup(selected_group_name)[ name_idx-1 ][ i ];
          std::vector<std::string>::iterator it = find( names.begin(),names.end(),elem );
          if( it != names.end() )
              param_lst.push_back( elem );
       }
   }
   // remove entries of given group from the entire list
   else if( name_idx < 0)
   {
       name_idx *= -1;
       param_lst = names;
       //for( int i=0; i<groups[selected_group_name][ name_idx-1].getLength();
       for( int i=0; i<groups.lookup(selected_group_name)[ name_idx-1].getLength();
               ++i )
       {
           vector<string>::iterator iter = find( param_lst.begin(), param_lst.end(),
                                                 //(const char *)groups[selected_group_name][ name_idx-1][ i ] );
                                                 (const char *)groups.lookup(selected_group_name)[ name_idx-1][ i ] );
           if (iter != param_lst.end())
               param_lst.erase( iter );
       }
   }
}



