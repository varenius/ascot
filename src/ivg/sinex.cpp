#include "sinex.h"
#include "param_list.h"

namespace ivg
{
    
    
// ...........................................................................
Sinex::Sinex()
// ...........................................................................
{
}
// ...........................................................................
Sinex::Sinex( const ivg::Sinex &other )
// ...........................................................................
{
    (*this) = other;
}
// ...........................................................................
Sinex & Sinex::operator=( const Sinex &other )
// ...........................................................................
{
#if DEBUG_VLBI >=2
   cerr << "+++ Sinex & Sinex::operator=( const Sinex &other )" << endl;
   tictoc tim;
   tim.tic();
#endif
        _path = other._path;
        
        _start = other._start;
        _end = other._end;
        
        _station_names = other._station_names;
        
        _src_assignment = other._src_assignment; // Pointer!!!
        _sta_assignment = other._sta_assignment;       
        
        _file_comment << other._file_comment.rdbuf();

        _parameter.resize( other._parameter.size() );
        copy( other._parameter.begin(), other._parameter.end(), _parameter.begin() );

        _discontinuities = other._discontinuities;
        _stats = other._stats;

        _N = other._N;
        _n_side = other._n_side;

#if DEBUG_VLBI >=2
   cerr << "--- Sinex & Sinex::operator=( const Sinex &other )" << " : " << tim.toc() << " s " << endl;
#endif 
   
    return *this;
}
// ...........................................................................   
Sinex::Sinex(string path ,string tropo_snx)
// ...........................................................................
{
#if DEBUG_VLBI >=2
    cerr<<"+++ Sinex::Sinex(string path, bool init_neq )" <<endl;
    tictoc tim;
    tim.tic();
#endif
    
    // full path to the sinex file
    _path = path;
   
    // data matrices
    ivg::Matrix N(1,1,0.0);
    ivg::Matrix n_side(1,1,0.0);
    vector<double> mjd_parameter;

    // intern save structure and relation between apriori and estimate block
    map<int,int> apriori_estimate_assignment;
    
    // need to know if parameter vector initialized or not
    // (because of unknown order of APRIORI and ESTIMATE block)
    bool parameter_vec_init=false;
    
    ifstream inStream(path.c_str(), ios::in);
    if( !inStream.is_open() ){
        throw runtime_error( "Sinex::Sinex(string path, bool init_neq ): Failed to open file: "+path );
    }else{

        string line;
        while (getline(inStream,line,'\n'))
        {
            if (line.substr(0,1)=="%"&&line.substr(1,7)!="=ENDTRO")
            {
                try{
                int y1 = atoi(line.substr(32,2).c_str());
                int y2 = atoi(line.substr(45,2).c_str());
                int doy1 = atoi(line.substr(35,3).c_str());
                int doy2 = atoi(line.substr(48,3).c_str());
                int sec1 = atoi(line.substr(39,5).c_str());
                int sec2 = atoi(line.substr(52,5).c_str());

                (y1<70)?y1 += 2000:y1 += 1900;
                (y2<70)?y2 += 2000:y2 += 1900;

                ivg::Date start(y1,(double) doy1+((double) sec1)/86400.0);
                ivg::Date end(y2,(double) doy2+((double) sec2)/86400.0);

                // setting start and epoch of data from sinex file
                _start = start;
                _end = end;

//                statistics NoPH; // Number of Parameters Header
//                NoPH.name = "NUMBER OF PARAMETERS HEADER";
//                NoPH.shorty = "NoPH";
//                NoPH.value = _stod(line.substr(60,5));
//                _stats[NoPH.shorty] = NoPH;

                }catch(exception dsrgs){
                    
                }
                

            }
            else if (line.find("+FILE/COMMENT")!=string::npos)
            {
                statistics tmp;

                while (getline(inStream,line,'\n')&&line.substr(0,1)!="-")
                {
                    if (line.find("WRMS")!=string::npos)
                    {
                        tmp.name = "WEIGHTED ROOT MEAN SQUARED";
                        tmp.shorty = "WRMS";
                        tmp.value = _stod(_DtoE(line.substr(32,line.size()-32)));
                        _stats[tmp.shorty] = tmp;
                    }
                    else if (line.find("RMS")!=string::npos)
                    {
                        tmp.name = "ROOT MEAN SQUARED";
                        tmp.shorty = "RMS";
                        tmp.value = _stod(_DtoE(line.substr(32,line.size()-32)));
                        _stats[tmp.shorty] = tmp;
                    }
                }
                // _file_comment <<line << endl;
                


              }else if(line.find("+TROP/STA_COORDINATES")!=string::npos){  
                while (getline(inStream, line, '\n') && line.substr(0,1) != "-"){
                    if (line.substr(0,1) != "*") {
                        string code = remove_spaces_end(line.substr(1,4));
                        boost::to_upper(code);
//                        string domes_no = remove_spaces_end(line.substr(9,9));
                        string ivs_name = remove_spaces_end(line.substr(1,4));
                        boost::to_upper(ivs_name);
//                        replace_string_in_place(ivs_name, " ", "_" );
                        map<ivg::staname,string> names;
                        names[ivg::staname::cdp] = code;
//                        names[ivg::staname::domes_no] = domes_no;
                        names[ivg::staname::ivs_name] = ivs_name;
                        
                        _sta_assignment[code] = names;
                        _sta_assignment[ivs_name] = names;
                    }
                }

            }else if(line.find("+SITE/ID")!=string::npos){  
                  
                while (getline(inStream, line, '\n') && line.substr(0,1) != "-"){
                    
                    if (line.substr(0,1) != "*") {
                        
                        string code = remove_spaces_end(line.substr(1,4));
                        string domes_no = remove_spaces_end(line.substr(9,9));
//                        string ivs_name = remove_spaces_end(line.substr(21,8));
//                        replace_string_in_place(ivs_name, " ", "_" );

                        map<ivg::staname,string> names;
                        names[ivg::staname::cdp] = code;
                        names[ivg::staname::domes_no] = domes_no;
//                        names[ivg::staname::ivs_name] = ivs_name;
                        
                        _sta_assignment[code] = names;
//                        _sta_assignment[ivs_name] = names;
                    }
                }
                
            }else if(line.find("+TROP/SOLUTION")!=string::npos || line.find("+SOLUTION/ESTIMATE")!=string::npos){

                // save which block we are parsing right now
                string solution;
//                if(line.find("+SOLUTION/APRIORI")!=string::npos)
//                    solution = "APRIORI";
//                else
                    solution = "ESTIMATE";
                int index = 0;

                
                bool paramAvailable = true;
                bool withSTD = false;
                std::size_t linepos;
                std::size_t lineposSTD;
                std::size_t lineposCode;
                std::size_t lineposTime;
                vector<string> columns;
                std::vector<string>::iterator  p;

                int zeile = 0;  //Nur erste Zeile nach "+TROP/SOLUTION" kann Spaltennamen enthalten, sonst Kommentar (z.B. "day49")
                
                while (getline(inStream,line,'\n')&&line.substr(0,1)!="-")
                {
                    zeile++;
                    
                    if (line.substr(0,1) == "*" && zeile==1) {

                            boost::char_separator<char> sep(" ");
                            boost::tokenizer<boost::char_separator<char>> tokens(line.substr(1, line.length() - 1), sep);
                            
                            for (const auto& t : tokens) {
                                columns.push_back(t);
                            }
                            
                            p = std::find(columns.begin(), columns.end(), tropo_snx);
                            if (p == columns.end() || (*p).empty()) {
//                                cerr << "not found------------" << endl;
                                paramAvailable = false;
                            } else {
//                                cerr << "found:" << *p << endl;
                                linepos = line.find(*p);
                                if (*(p + 1) == "STDEV") {
                                    withSTD = true;
                                    lineposSTD = line.find(*(p + 1));
                                }
                                lineposCode =  line.find("SITE");
                                lineposTime =  line.find(columns.at(1));
                            }
//                            cerr<<"endheader"<<endl;
                        
                    }
                    else if (line.substr(0,1)!="*" && paramAvailable)
                    {
                        string code = line.substr(lineposCode,4).c_str();
                        boost::to_upper(code);

                        ivg::Date epoch(line.substr(lineposTime,columns.at(1).size()),"SINEX");

                        double TROTOT, TROTOT_std;
                        stringstream sstr1;
                        sstr1 << remove_spaces_end(_DtoE(line.substr(linepos, 6)));
                        if(withSTD){
                            stringstream sstr2;
                            sstr2 << remove_spaces_end(_DtoE(line.substr(lineposSTD, 5)));
                            sstr2 >> TROTOT_std;
                        }
                        
                        // check if parameter vector already has been initialized or not
                        // depends on if the first block is ESTIMATE or APRIORI
                        if(!parameter_vec_init)
                        {
//                            vector<string>::const_iterator pos_iter = find( paramtype_snx.begin(), paramtype_snx.end(), tropo_snx );
                            vector<string>::const_iterator pos_iter = find( paramtype_snx.begin(), paramtype_snx.end(), "TROTOT" );
                            int pos = pos_iter - paramtype_snx.begin();
                                                        
                            // we use a struct to save parameter specific information just temporary
                            struct param_init{
                                ivg::paramtype type;
                                string name;
                                int order;
                            } p_tmp;
                            
                            // CRITICAL TO USE HARD CODED NUMBERS HERE !!!!!!!!
                            //   0       1       2       3       4         5         6        7      8      9      10        11       12      13        14       15        16
                           //{ "STAX", "STAY", "STAZ", "CLO", "TROTOT", "TGNTOT", "TGETOT", "XPO", "YPO", "UT1", "NUT_X", "NUT_Y", "RS_RA", "RS_DE", "CL_BR", "NUT_LN", "NUT_OB" } );
                            p_tmp = {(ivg::paramtype)pos, _sta_assignment[code][ivg::staname::ivs_name], 0 };
                            
                            ivg::Param tmp_param = ivg::Param(p_tmp.type, p_tmp.name, epoch, TROTOT, p_tmp.order);
                            if(tropo_snx == "TROTOT"){
                                tmp_param.set_apriori(TROTOT);
                            }else if(tropo_snx == "TROWET"){
                                tmp_param.set_estimate(TROTOT);
                            }
                            if(withSTD){
                                tmp_param.set_standard_deviation(TROTOT_std);
                            }
                            _parameter.push_back(tmp_param);

                            apriori_estimate_assignment[index] = _parameter.size()-1;
                        }
                        // if parameter vector already has been initialized
                        else
                        {
                            int vec_position = apriori_estimate_assignment[index];
                            
                            if(tropo_snx == "TROTOT"){
                                _parameter.at(vec_position).set_apriori(TROTOT);
                            }else if(tropo_snx == "TROWET"){
                                _parameter.at(vec_position).set_estimate(TROTOT);
                            }
                            if(withSTD){
                                _parameter.at(vec_position).set_standard_deviation(TROTOT_std);
                            }
                        }
                        index+=1;
                    }
                }

                    if (paramAvailable) 
                    {
                        statistics tmp; // Number of Parameters Apriori Block
                        tmp.name = "NUMBER OF PARAMETERS APRIBLOCK";
                        tmp.shorty = "NoPAB";
                        tmp.value = _parameter.size() / 1.0;
                        _stats[tmp.shorty] = tmp;

//                        _parameter.data()->show();

                        parameter_vec_init = true;
                    }
                }
                else if(line.substr(0,line.size()) == "+SOLUTION/CONSTRAINT_INFO"){
                
                double nCons=0;
                
                while (getline(inStream, line, '\n') && line.substr(0,1) != "-"){
                    if (line.substr(0,1) != "*") 
                       nCons++;
                }

                statistics tmp;
                tmp.name = "NUMBER OF CONSTRAINTS";
                tmp.shorty = "NoC";
                tmp.value = nCons;
                _stats[tmp.shorty] = tmp;

                _stats["NoO+NoC"].value = nCons+_stats["NoO"].value;

            }
        }
    }
    
    inStream.close();
    
#if DEBUG_VLBI >=2
    cerr<<"--- Sinex::Sinex(string path, bool init_neq )"
            <<" : "<<tim.toc()<<" s "<<endl;
#endif  
}
// .......................
// ...........................................................................   
Sinex::Sinex(string path, bool init_neq )
// ...........................................................................
{
#if DEBUG_VLBI >=2
    cerr<<"+++ Sinex::Sinex(string path, bool init_neq )" <<endl;
    tictoc tim;
    tim.tic();
#endif
        
    // full path to the sinex file
    _path = path;
   
    // data matrices
    ivg::Matrix N(1,1,0.0);
    ivg::Matrix n_side(1,1,0.0);
    vector<double> mjd_parameter;

    // intern save structure and relation between apriori and estimate block
    map<int,int> apriori_estimate_assignment;
    
    // need to know if parameter vector initialized or not
    // (because of unknown order of APRIORI and ESTIMATE block)
    bool parameter_vec_init=false;
    
    ifstream inStream(path.c_str(), ios::in);
    if( !inStream.is_open() ){
        throw runtime_error( "Sinex::Sinex(string path, bool init_neq ): Failed to open file: "+path );
    }else{

        string line;
        while (getline(inStream,line,'\n'))
        {
//            cerr << line << endl;
            if (line.substr(0,1)=="%"&&line.substr(1,6)!="ENDSNX")
            {

                int y1 = atoi(line.substr(32,2).c_str());
                int y2 = atoi(line.substr(45,2).c_str());
                int doy1 = atoi(line.substr(35,3).c_str());
                int doy2 = atoi(line.substr(48,3).c_str());
                int sec1 = atoi(line.substr(39,5).c_str());
                int sec2 = atoi(line.substr(52,5).c_str());

                (y1<70)?y1 += 2000:y1 += 1900;
                (y2<70)?y2 += 2000:y2 += 1900;

                ivg::Date start(y1,(double) doy1+((double) sec1)/86400.0);
                ivg::Date end(y2,(double) doy2+((double) sec2)/86400.0);

                // setting start and epoch of data from sinex file
                _start = start;
                _end = end;

                statistics NoPH; // Number of Parameters Header
                NoPH.name = "NUMBER OF PARAMETERS HEADER";
                NoPH.shorty = "NoPH";
                NoPH.value = _stod(line.substr(60,5));
                _stats[NoPH.shorty] = NoPH;

            }
            else if (line.find("+FILE/COMMENT")!=string::npos)
            {
                statistics tmp;
                
                while (getline(inStream,line,'\n')&&line.substr(0,1)!="-")
                {
                    if (line.find("WRMS")!=string::npos)
                    {
                        tmp.name = "WEIGHTED ROOT MEAN SQUARED";
                        tmp.shorty = "WRMS";
                        tmp.value = _stod(_DtoE(line.substr(32,line.size()-32)));
                        _stats[tmp.shorty] = tmp;
                    }
                    else if (line.find("RMS")!=string::npos)
                    {
                        tmp.name = "ROOT MEAN SQUARED";
                        tmp.shorty = "RMS";
                        tmp.value = _stod(_DtoE(line.substr(32,line.size()-32)));
                        _stats[tmp.shorty] = tmp;
                    }
                }
                // _file_comment <<line << endl;
                
            }
            else if (line.find("+SOLUTION/STATISTICS")!=string::npos)
            {
                statistics tmp;

                while (getline(inStream,line,'\n')&&line.substr(0,1)!="-")
                {

                    if (line.find("NUMBER OF DEGREES OF FREEDOM")!=string::npos)
                    { //NoDoF
                        tmp.name = "NUMBER OF DEGREES OF FREEDOM";
                        tmp.shorty = "NoDoF";
                        tmp.value = _stod(_DtoE(line.substr(32,line.size()-32)));
                        _stats[tmp.shorty] = tmp;
                    }
                    else if (line.find("SQUARE SUM OF O-C")!=string::npos)
                    { //LTPL
                        tmp.name = "SQUARE SUM OF O-C";
                        tmp.shorty = "LTPL";
                        tmp.value = _stod(_DtoE(line.substr(32,line.size()-32)));
                        _stats[tmp.shorty] = tmp;
                    }
                    else if (line.find("NUMBER OF OBSERVATIONS")!=string::npos)
                    { //NoO
                        tmp.name = "NUMBER OF OBSERVATIONS";
                        tmp.shorty = "NoO";
                        tmp.value = _stod(_DtoE(line.substr(32,line.size()-32)));
                        _stats[tmp.shorty] = tmp;

                        tmp.name = "NUMBER OF OBSERVATIONS+CONSTS";
                        tmp.shorty = "NoO+NoC";
                        tmp.value = _stod(_DtoE(line.substr(32,line.size()-32)));
                        _stats[tmp.shorty] = tmp;
                    }
                    else if (line.find("WRMS OF POSTFIT RESIDUALS")!=string::npos)
                    { //WoPR
                        tmp.name = "WRMS OF POSTFIT RESIDUALS";
                        tmp.shorty = "WoPR";
                        tmp.value = _stod(_DtoE(line.substr(32,line.size()-32)));
                        _stats[tmp.shorty] = tmp;
                    }
                    else if (line.find("NUMBER OF UNKNOWNS")!=string::npos)
                    { //NoU
                        tmp.name = "NUMBER OF UNKNOWNS";
                        tmp.shorty = "NoU";
                        tmp.value = _stod(_DtoE(line.substr(32,line.size()-32)));
                        _stats[tmp.shorty] = tmp;
                    }
                    else if (line.find("VARIANCE FACTOR")!=string::npos)
                    { //VF
                        tmp.name = "VARIANCE FACTOR";
                        tmp.shorty = "VF";
                        tmp.value = _stod(_DtoE(line.substr(32,line.size()-32)));
                        _stats[tmp.shorty] = tmp;
                    }
                    else if (line.find("SQUARE SUM OF RESIDUALS")!=string::npos)
                    { //VTPV
                        tmp.name = "SQUARE SUM OF RESIDUALS";
                        tmp.shorty = "VTPV";
                        tmp.value = _stod(_DtoE(line.substr(32,line.size()-32)));
                        _stats[tmp.shorty] = tmp;
                    }

                }

              }else if(line.find("+SITE/ID")!=string::npos){  
                  
                while (getline(inStream, line, '\n') && line.substr(0,1) != "-"){
                    
                    if (line.substr(0,1) != "*") {
                        
                        string code = remove_spaces_end(line.substr(1,4));
                        string domes_no = remove_spaces_end(line.substr(9,9));
                        string ivs_name = remove_spaces_end(line.substr(21,8));
                        replace_string_in_place(ivs_name, " ", "_" );

                        map<ivg::staname,string> names;
                        names[ivg::staname::cdp] = code;
                        names[ivg::staname::domes_no] = domes_no;
                        names[ivg::staname::ivs_name] = ivs_name;
                        
                        _sta_assignment[code] = names;
                        _sta_assignment[ivs_name] = names;
                        
                        _station_names.push_back(ivs_name);
                    }
                }
                
            }else if(line.find("+SOLUTION/EPOCHS")!=string::npos){ 
                
                while (getline(inStream, line, '\n') && line.substr(0,1) != "-"){
                    if (line.substr(0,1) != "*") {
                        
                        string code = remove_spaces_end(line.substr(1,4));
                        int soln =  stoi(remove_spaces_end(line.substr(9,4)));
                        
                        ivg::Date epoch_start( remove_spaces_end(line.substr(16,12)), "SINEX" );
                        ivg::Date epoch_end( remove_spaces_end(line.substr(29,12)), "SINEX" );
                        ivg::Date epoch_mean( remove_spaces_end(line.substr(42,12)), "SINEX" );
                        
                        
                        _discontinuities[code].push_back(epoch_start);
                    }
                }
                
            }else if(line.find("+SOURCE/ID")!=string::npos){ 
                
                while (getline(inStream, line, '\n') && line.substr(0,1) != "-"){
                    if (line.substr(0,1) != "*") {
                        
                        //assignment between code and iers name of source
                        string iers_name = remove_spaces_end(line.substr(6,8));
                        string ivs_name;
                        // in case of gfz-snx files we don't have IAU names
                        if( path.find("gfz") != string::npos )
                            ivs_name = remove_spaces_end(line.substr(32,8));  
                        // in case of dgf files, the ivs name is located somewhere else
                        else if( path.find("dgf") != string::npos )
                              ivs_name = remove_spaces_end(line.substr(34));  
                        // all other cases
                        else 
                            ivs_name = remove_spaces_end(line.substr(43,8));

                        map<ivg::srcname,string> names;
                        names[ivg::srcname::iers] = iers_name;
                        names[ivg::srcname::ivs] = ivs_name;
                        
                        _src_assignment[line.substr(1,4)] = names;
                        _src_assignment[ivs_name] = names;

                    }
                }
                
            }else if(line.find("+SOLUTION/APRIORI")!=string::npos || line.find("+SOLUTION/ESTIMATE")!=string::npos){

                // save which block we are parsing right now
                string solution;
                if(line.find("+SOLUTION/APRIORI")!=string::npos)
                    solution = "APRIORI";
                else
                    solution = "ESTIMATE";
                
                while (getline(inStream,line,'\n')&&line.substr(0,1)!="-")
                {
                    if (line.substr(0,1)!="*")
                    {
                        // get information from file
                        int index = stoi(line.substr(1,5));
                        string type = remove_spaces_end(line.substr(7,6));
                        string code = line.substr(14,4).c_str();
                        // SBIN can be used as stacked
                        int stacked = stoi(line.substr(22,4));
                        ivg::Date epoch(line.substr(27,12),"SINEX");
                        mjd_parameter.push_back(epoch.get_double_mjd());
                        
                        double apri_or_esti, const_or_std;
                        
                        stringstream sstr1; 
                        sstr1 << remove_spaces_end(_DtoE(line.substr(47,21))); 
                        sstr1 >> apri_or_esti;
                        
                        stringstream sstr2;
                        sstr2 << remove_spaces_end(_DtoE(line.substr(69,11))); 
                        sstr2 >> const_or_std;
                        
                        // check if parameter vector already has been initialized or not
                        // depends on if the first block is ESTIMATE or APRIORI
                        if(!parameter_vec_init)
                        {
                            vector<string>::const_iterator pos_iter = find( paramtype_snx.begin(), paramtype_snx.end(), type );
                            int pos = pos_iter - paramtype_snx.begin();
                            
                            // we use a struct to save parameter specific information just temporary
                            struct param_init{
                                ivg::paramtype type;
                                string name;
                                int order;
                                int stacked;
                            } p_tmp;
                            
                            // CRITICAL TO USE HARD CODED NUMBERS HERE !!!!!!!!
                            //   0       1       2       3       4         5         6        7      8      9      10        11       12      13        14       15        16
                           //{ "STAX", "STAY", "STAZ", "CLO", "TROTOT", "TGNTOT", "TGETOT", "XPO", "YPO", "UT1", "NUT_X", "NUT_Y", "RS_RA", "RS_DE", "CL_BR", "NUT_LN", "NUT_OB" } );
                            if(pos >= 0 && pos <=2 )
                            {
                                p_tmp = {(ivg::paramtype)pos, _sta_assignment[code][ivg::staname::cdp], 0 };
                            }
                            else if(pos >=4 && pos <= 6)
                                p_tmp = {(ivg::paramtype)pos, _sta_assignment[code][ivg::staname::cdp], 0 };
                            else if(pos >=7 && pos <=11)
                                p_tmp = {(ivg::paramtype)pos, "EOP", 0 };
                            else if(pos >= 12 && pos <=13 )
                                p_tmp = {(ivg::paramtype)pos, _src_assignment[code][ivg::srcname::ivs], 0 };
                            else if(pos >= 14)
                            {
                                if(type == "VELX" || type == "VELY" || type == "VELZ" )
                                {
                                    
                                    if(type == "VELX")
                                        p_tmp = {ivg::paramtype::stax, _sta_assignment[code][ivg::staname::cdp], 1 };
                                    else if(type == "VELY")
                                        p_tmp = {ivg::paramtype::stay, _sta_assignment[code][ivg::staname::cdp], 1 };
                                    else if(type == "VELZ")
                                        p_tmp = {ivg::paramtype::staz, _sta_assignment[code][ivg::staname::cdp], 1 }; 
                                   
                                }
                                else if(type == "XPOR")
                                   p_tmp = {ivg::paramtype::xpo, "EOP", 1 };
                                else if(type == "YPOR")
                                   p_tmp = {ivg::paramtype::ypo, "EOP", 1 };
                                else if(type == "UT")
                                   p_tmp = {ivg::paramtype::ut1, "EOP", 0 };
                                else if(type == "LOD")
                                   p_tmp = {ivg::paramtype::ut1, "EOP", 1 };
                                else if(type == "NUT_LN")
                                     p_tmp = {ivg::paramtype::nutln, "EOP", 0 };
                                else if(type == "NUT_OB")
                                      p_tmp = {ivg::paramtype::nutob, "EOP", 0 };
                                else if(type.find("CLO") != string::npos)
                                    p_tmp = {ivg::paramtype::clo, _sta_assignment[code][ivg::staname::cdp], stoi(type.substr(3)) };
                                else if(type.find("CL_BR") != string::npos)
                                    p_tmp = {ivg::paramtype::cbr, "----", 0 };                                
                                else
                                    throw runtime_error("void _read_snx(ivg::Session *, Setting *, const string ): Unknown paramtype "+type+" in "+path);
                            }
                            
                            // adjust leap seconds in case of UT / UT1
//                            if(p_tmp.type == ivg::paramtype::ut1 && p_tmp.order == 0)
//                            {
//                                // if "UT" is smaller than 1500.0, we need to add leap seconds to be consistent 
//                                if(apri_or_esti < 0 && apri_or_esti > -1500.0)
//                                {
//                                    log<WARNING>("!!! Correcting for leap seconds in _param_list for "+path);
//                                    apri_or_esti -= epoch.get_leap_sec()*1000;
//                                }
//                                else if(apri_or_esti >= 0 && apri_or_esti < 1500.0)
//                                {
//                                    log<WARNING>("!!! Correcting for leap seconds in _param_list for "+path);
//                                    apri_or_esti += epoch.get_leap_sec()*1000;
//                                }
//                            }
                            
                            // dgf uses always the same epoch (1.1.2000) for a source in every available snx file
                            // in order to be able to set up sources as local parameter, we need to assign a daily epoch
//                            if(path.find("dgf") != string::npos && (p_tmp.type == ivg::paramtype::ra || p_tmp.type == ivg::paramtype::dec) )
//                            {
//                                double mid_session_mjd = (_start.get_double_mjd() + _end.get_double_mjd() )/2.0;
//                                epoch = ivg::Date(mid_session_mjd);
//                            }
                            
                            
                            if(solution == "APRIORI")
                            {
                                _parameter.push_back(ivg::Param(p_tmp.type, p_tmp.name, epoch, apri_or_esti, p_tmp.order));
                                _parameter.back().set_stacked(p_tmp.stacked);
                            }
                            else if(solution == "ESTIMATE")
                            {
                                ivg::Param tmp_param = ivg::Param(p_tmp.type, p_tmp.name, epoch, 0.0, p_tmp.order);
                                tmp_param.set_estimate(apri_or_esti);
                                tmp_param.set_standard_deviation(const_or_std);
                                tmp_param.set_stacked(p_tmp.stacked);
                                _parameter.push_back(tmp_param);
                            }

                            apriori_estimate_assignment[index] = _parameter.size()-1;
                        }
                        // if parameter vector already has been initialized
                        else
                        {
                            int vec_position = apriori_estimate_assignment[index];
                            
                            if(solution == "APRIORI")
                            {
                                _parameter.at(vec_position).set_apriori(apri_or_esti);
                                _parameter.at(vec_position).set_stacked(stacked);
                            }
                            else if(solution == "ESTIMATE")
                            {
                                _parameter.at(vec_position).set_estimate(apri_or_esti);
                                _parameter.at(vec_position).set_standard_deviation(const_or_std);
                                _parameter.at(vec_position).set_stacked(stacked);
                            }
                        }
                    }
                }
                
            statistics tmp; // Number of Parameters Apriori Block
            tmp.name = "NUMBER OF PARAMETERS APRIBLOCK";
            tmp.shorty = "NoPAB";
            tmp.value = _parameter.size()/1.0;
            _stats[tmp.shorty] = tmp;

            parameter_vec_init = true;
            }
                else if(line.substr(0,line.size()) == "+SOLUTION/CONSTRAINT_INFO"){
                
                double nCons=0;
                
                while (getline(inStream, line, '\n') && line.substr(0,1) != "-"){
                    if (line.substr(0,1) != "*") 
                       nCons++;
                }

                statistics tmp;
                tmp.name = "NUMBER OF CONSTRAINTS";
                tmp.shorty = "NoC";
                tmp.value = nCons;
                _stats[tmp.shorty] = tmp;

                _stats["NoO+NoC"].value = nCons+_stats["NoO"].value;

            }
            else if (line.find("+SOLUTION/NORMAL_EQUATION_MATRIX")!=string::npos||line.find("+SOLUTION/DECOMPOSED_NORMAL_MATRIX")!=string::npos)
            {
                // if we dont want the matrices, we stop reading the file
                if(init_neq == false)
                {
                    inStream.close();
                    break;
                }

                _N.resize(_parameter.size(),_parameter.size(),0.0);

                while (getline(inStream,line,'\n')&&line.substr(0,1)!="-")
                {
                    if (line.substr(0,1)!="*")
                    {

                        int n = stoi(line.substr(1,5));
                        int m = stoi(line.substr(7,5));
                        int lineSize = line.size();
                        
                        if (lineSize==34 || lineSize==35)
                        {
                            _N(n-1,m-1) = _stod(_DtoE(line.substr(13,21)));
                            _N(m-1,n-1) = _N(n-1,m-1);
                        }
                        else if (lineSize==56)
                        {
                            _N(n-1,m-1) = _stod(_DtoE(line.substr(13,21)));
                            _N(n-1,m) = _stod(_DtoE(line.substr(35,21)));
                            _N(m-1,n-1) = _N(n-1,m-1);
                            _N(m,n-1) = _N(n-1,m);
                        }
                        else if (lineSize==78)
                        {

                            _N(n-1,m-1) = _stod(_DtoE(line.substr(13,21)));
                            _N(n-1,m) = _stod(_DtoE(line.substr(35,21)));
                            _N(n-1,m+1) = _stod(_DtoE(line.substr(57,21)));

                            _N(m-1,n-1) = _N(n-1,m-1);
                            _N(m,n-1) = _N(n-1,m);
                            _N(m+1,n-1) = _N(n-1,m+1);
                        }
                        else
                        {
                            throw runtime_error("void _read_snx(ivg::Session *, Setting *, const string ): Unexpected line length of matrix block in "+path);
                        }
                    }
                }

            }
            else if (line.find("+SOLUTION/NORMAL_EQUATION_VECTOR")!=string::npos||line.find("+SOLUTION/DECOMPOSED_NORMAL_VECTOR")!=string::npos)
            {
                _n_side.resize(_parameter.size(),1,0.0);

                while (getline(inStream, line, '\n') && line.substr(0,1) != "-"){
                    if (line.substr(0,1) != "*"){
                        int n = stoi(line.substr(1,5));
                        
                        // in some cases there is a different type of VECTOR-Block
                        if(line.size() > 50)
                            _n_side(n-1,0) = _stod(_DtoE(line.substr(47,23)));
                        else
                            _n_side(n-1,0) = _stod(_DtoE(line.substr(7,23)));
                    }
                }
            }
        }
    }
        
    inStream.close();
    
#if DEBUG_VLBI >=2
    cerr<<"--- Sinex::Sinex(string path, bool init_neq )"
            <<" : "<<tim.toc()<<" s "<<endl;
#endif  
}
// ...........................................................................
ivg::Trf Sinex::get_trf(reftype type)
// ...........................................................................
{
#if DEBUG_VLBI >=2
    cerr<<"+++ ivg::Trf Sinex::get_trf(reftype type)" <<endl;
    tictoc tim;
    tim.tic();
#endif
    
    // vector of analysis_stations
    vector<ivg::Analysis_station> stations;
       
    // create trf based on stations vector
    double mid_session_mjd = (_start.get_double_mjd() + _end.get_double_mjd() )/2.0;
    ivg::Trf trf = ivg::Trf(_path, ivg::Date(mid_session_mjd), stations);
    
    struct disconti{
        ivg::Date refepoch;
        ivg::Matrix pos0;
        ivg::Matrix vel0;
        vector<ivg::Date> epochs;
    };
    
    map< string, int > soln_counter;
    
    // go through ALL parameter
    for(vector<ivg::Param>::iterator param_iter = _parameter.begin(); param_iter!= _parameter.end(); ++param_iter)
    {
        // only if the paramter is a X-coordinate
        if(param_iter->get_type() == ivg::paramtype::stax && param_iter->get_order() == 0)
        {
            //only if X,Y,Z are in ascending parameter-order
            if((param_iter+1)->get_type() == ivg::paramtype::stay && (param_iter+2)->get_type() == ivg::paramtype::staz )
            {
                // only if epochs of X,Y,Z are equal
                if( (param_iter)->get_epoch() == (param_iter+1)->get_epoch() && (param_iter+1)->get_epoch() == (param_iter+2)->get_epoch() )
                {   
                    ivg::Matrix xyz0; 
                    ivg::Matrix xyz0_std(3,1,0.0);
                    ivg::Matrix vel0(3,1,0.0);
                    ivg::Matrix vel0_std(3,1,0.0);
                    if( type == reftype::estimate )
                    {
                        xyz0 = ivg::Matrix(vector<double>{(param_iter)->get_estimate(),(param_iter+1)->get_estimate(),(param_iter+2)->get_estimate()});
                        xyz0_std = ivg::Matrix(vector<double>{(param_iter)->get_standard_deviation(),
                                                              (param_iter+1)->get_standard_deviation(),
                                                              (param_iter+2)->get_standard_deviation()});
                    }
                    else if( type == reftype::apriori )
                        xyz0 = ivg::Matrix(vector<double>{(param_iter)->get_apriori(),(param_iter+1)->get_apriori(),(param_iter+2)->get_apriori()});
                        
                    // we need to check if velocities are existent in _param_list
                    if(  (param_iter+3) != _parameter.end() && param_iter->get_name() == (param_iter+3)->get_name() && (param_iter+3)->get_type() == ivg::paramtype::stax && (param_iter+3)->get_order() == 1 )
                    {
                        if( type == reftype::estimate )
                        {
                            vel0 = ivg::Matrix(vector<double>{(param_iter+3)->get_estimate(),(param_iter+4)->get_estimate(),(param_iter+5)->get_estimate()});
                            vel0_std = ivg::Matrix(vector<double>{(param_iter+3)->get_standard_deviation(),
                                                                  (param_iter+4)->get_standard_deviation(),
                                                                  (param_iter+5)->get_standard_deviation()});
                        }                   
                        else if( type == reftype::apriori )
                            vel0 = ivg::Matrix(vector<double>{(param_iter+3)->get_apriori(),(param_iter+4)->get_apriori(),(param_iter+5)->get_apriori()});
                    }

                    ivg::Date epoch = param_iter->get_epoch();
                    string cdp = param_iter->get_name();
                    
                    ivg::Analysis_station * sta_iter;
                    // we need the corresponding analysis_station to the param_iter
                    if(trf.get_station(&sta_iter,param_iter->get_name(), ivg::staname::description))
                    {
                        string code = _sta_assignment[cdp][staname::cdp];
                        
                        sta_iter->add_discontinuity( xyz0, vel0, epoch, 
                                                     _discontinuities[code].at(soln_counter[code]),
                                                     xyz0_std, vel0_std );
                        soln_counter[code]++;
                    }
                    else
                    {
                        vector<ivg::Date> single_discont = {ivg::Date(1970,1.0)};
                        ivg::Analysis_station sta( xyz0, vel0, xyz0_std, vel0_std, epoch, single_discont , _sta_assignment[cdp]);
                        trf.push_back(sta);
                        
                        soln_counter[_sta_assignment[cdp][staname::cdp]]++;
                    }
                }
                else
                    throw runtime_error( "ivg::Trf Sinex::get_trf(reftype type): Epochs of X,Y,Z of station "+param_iter->get_name()+" not equal.");
            }
            else
                throw runtime_error( "ivg::Trf Sinex::get_trf(reftype type): X,Y,Z of station "+param_iter->get_name()+" not in correct order.");
        }
    }
 
    
#if DEBUG_VLBI >=2
    cerr<<"--- ivg::Trf Sinex::get_trf(reftype type)"
            <<" : "<<tim.toc()<<" s "<<endl;
#endif  
    
    return trf;
}
// ...........................................................................
ivg::Crf Sinex::get_crf(reftype type)
// ...........................................................................
{
#if DEBUG_VLBI >=2
    cerr<<"+++ ivg::Crf Sinex::get_crf(reftype type)" <<endl;
    tictoc tim;
    tim.tic();
#endif
    
    // vector of analysis_stations
    vector<ivg::Source> sources;
       
    // create crf based on source vector (empty)
    ivg::Crf crf = ivg::Crf(_path, sources );
    
    // go through ALL parameter
    for(vector<ivg::Param>::iterator param_iter = _parameter.begin(); param_iter!= _parameter.end(); ++param_iter)
    {
        // only if the paramter is right ascension
        if(param_iter->get_type() == ivg::paramtype::ra && param_iter->get_order() == 0)
        {
            //only if declination is in ascending parameter-order
            if((param_iter+1)->get_type() == ivg::paramtype::dec  )
            {
                // only if epochs of ra and dec are equal
                if( (param_iter)->get_epoch() == (param_iter+1)->get_epoch() )
                {   
                    double ra0, dec0;
                    if( type == reftype::estimate )
                    {
                        ra0 = (param_iter)->get_estimate();
                        dec0 = (param_iter+1)->get_estimate();
                    }
                    else if( type == reftype::apriori )
                    {
                        ra0 = (param_iter)->get_apriori();
                        dec0 = (param_iter+1)->get_apriori();
                    }

                    ivg::Date epoch = param_iter->get_epoch();
                    string ivs_name = param_iter->get_name();

                    ivg::Source * src_iter;
                    // we need the corresponding source to the param_iter
                    // if source is already existing
                    if(crf.get_source(&src_iter,param_iter->get_name()))
                    {
                        src_iter->add_local_position(ra0, dec0, epoch);
                    }
                    else
                    {
                        ivg::Source src(ivg::srcname::ivs, ivs_name, ra0, dec0);
                        src.set_refepoch( epoch );
                        src.set_name(ivg::srcname::iers, _src_assignment[ivs_name][ivg::srcname::iers]);
                        crf.push_back(src);
                    }
                }
                else
                    throw runtime_error( "ivg::Crf Sinex::get_crf(reftype type): Epochs of ra and dec of source "+param_iter->get_name()+" not equal.");
            }
            else
                throw runtime_error( "ivg::Crf Sinex::get_crf(reftype type): ra and dec of source "+param_iter->get_name()+" not in correct order.");
        }
    }
    
#if DEBUG_VLBI >=2
    cerr<<"--- ivg::Crf Sinex::get_crf(reftype type)"
            <<" : "<<tim.toc()<<" s "<<endl;
#endif  
    
    return crf;
}
// ...........................................................................
ivg::Eop_series Sinex::get_eop_series(reftype type)
// ...........................................................................
{
#if DEBUG_VLBI >=2
    cerr<<"+++ ivg::Eop_series Sinex::get_eop_series(reftype type)" <<endl;
    tictoc tim;
    tim.tic();
#endif
        
    // data storing structure
    map< double, ivg::Matrix > data_series;
    
    double last_erp_mjd = 0.0;
    // go through ALL parameter
    for(vector<ivg::Param>::iterator param_iter = _parameter.begin(); param_iter!= _parameter.end(); ++param_iter)
    {
        if(param_iter->is_type({xpo, ypo, ut1, nutx, nuty},{0, 1}))
        {
            double mjd = param_iter->get_epoch().get_double_mjd();
            
            // save the last XPO, YPO, or UT1 epoch for the NUT epoch
            // NUT epoch is set to XPO,YPO,UT1 epoch
            if(param_iter->is_type({xpo, ypo, ut1},{0, 1}))
                last_erp_mjd = mjd;
            // if NUT is the first parameter we use this epoch
            else if(last_erp_mjd == 0.0 && param_iter->is_type({nutx, nuty},{0}))
                mjd = mjd;
            else
                mjd = last_erp_mjd;
                    
            if(data_series[mjd].cols() == 0)
            {
                data_series[mjd] = ivg::Matrix(1,17,0.0);
                data_series[mjd](0,0) = mjd;
            }
            
            double value, std;
            if( type == reftype::estimate )
                value = (param_iter)->get_estimate();
            else if( type == reftype::apriori )
                value = (param_iter)->get_apriori();
            
            std = (param_iter)->get_standard_deviation();
            
            // XPO
            if(param_iter->get_type() == ivg::paramtype::xpo && param_iter->get_order() == 0)
            {
                data_series[mjd](0,1) = value * ivg::mas2rad;
                data_series[mjd](0,9) = std * ivg::mas2rad;
            }
            // YPO
            else if(param_iter->get_type() == ivg::paramtype::ypo && param_iter->get_order() == 0)
            {
                data_series[mjd](0,2) = value * ivg::mas2rad;
                data_series[mjd](0,10) = std * ivg::mas2rad;
            }
            // UT1
            else if(param_iter->get_type() == ivg::paramtype::ut1 && param_iter->get_order() == 0)
            {
                data_series[mjd](0,3) = value * 1e-3 * ivg::s2rad;
                data_series[mjd](0,11) = std * 1e-3 * ivg::s2rad;
            }
            // XPO-RATE
            if(param_iter->get_type() == ivg::paramtype::xpo && param_iter->get_order() == 1)
            {
                data_series[mjd](0,4) = value * ivg::mas2rad;
                data_series[mjd](0,12) = std * ivg::mas2rad;
            }
            // YPO-RATE
            else if(param_iter->get_type() == ivg::paramtype::ypo && param_iter->get_order() == 1)
            {
                data_series[mjd](0,5) = value * ivg::mas2rad;
                data_series[mjd](0,13) = std * ivg::mas2rad;
            }
            // UT1-RATE
            else if(param_iter->get_type() == ivg::paramtype::ut1 && param_iter->get_order() == 1)
            {
                data_series[mjd](0,6) = value * 1e-3 * ivg::s2rad;
                data_series[mjd](0,14) = std * 1e-3 * ivg::s2rad;
            }
            // NUT_X
            else if(param_iter->get_type() == ivg::paramtype::nutx && param_iter->get_order() == 0)
            {
                data_series[mjd](0,7) = value * ivg::mas2rad;
                data_series[mjd](0,15) = std * ivg::mas2rad;
            }
            // NUT_Y
            else if(param_iter->get_type() == ivg::paramtype::nuty && param_iter->get_order() == 0)
            {
                data_series[mjd](0,8) = value * ivg::mas2rad;
                data_series[mjd](0,16) = std * ivg::mas2rad;
            }
        }
    }
    
    ivg::Matrix data;
    for(auto &data_row: data_series)
        data.append_rows(data_row.second);
    
#if DEBUG_VLBI >=2
    cerr<<"--- ivg::Eop_series Sinex::get_eop_series(reftype type)"
            <<" : "<<tim.toc()<<" s "<<endl;
#endif  
    
    // in case of unexpected eop data, return default eop series
    if(data.cols() != 17 || data(0,0) == 0.0)
        return ivg::Eop_series();
    else
        return ivg::Eop_series(_path, data);
}
// ...........................................................................
double Sinex::get_statistics(  std::string shorty )
// ...........................................................................
{
#if DEBUG_VLBI >=2
    cerr<<"+++ ivg::Eop_series Sinex::get_statistics(  std::string shorty )" <<endl;
    tictoc tim;
    tim.tic();
#endif

    return _stats[ shorty ].value;
    
#if DEBUG_VLBI >=2
    cerr<<"--- ivg::Eop_series Sinex::get_statistics(  std::string shorty )"
            <<" : "<<tim.toc()<<" s "<<endl;
#endif  
}

// ...........................................................................
map< string, ivg::Matrix > Sinex::get_tropo_delays( reftype type )
// ...........................................................................
{
#if DEBUG_VLBI >=2
    cerr << "+++ ivg::Trf Sinex::get_tropo_delays( reftype type )" << endl;
    tictoc tim;
    tim.tic();
#endif
    
    map< string, ivg::Matrix > tropo;
    ivg::Param_list param_list( _parameter );
        
    if( type == reftype::estimate )
        tropo = param_list.get_station_dependent_param( ivg::paramtype::zwd, ivg::param_def::totals, false, true );
    else if( type == reftype::apriori )
        tropo = param_list.get_station_dependent_param( ivg::paramtype::zwd, ivg::param_def::aprioris, false, true );
    else if( type == reftype::delta )
        tropo = param_list.get_station_dependent_param( ivg::paramtype::zwd, ivg::param_def::estimates, false, true );
    
    // change station names from CDP (which is used in SINEX) to IVS names (standard in ivg::ASCOT)    
    for( auto& it : tropo )
    {
        ivg::Matrix tmp = tropo[ it.first ];
        tropo[ _sta_assignment[it.first][ivg::staname::ivs_name] ] = tmp;
        tropo.erase( _sta_assignment[it.first][ivg::staname::cdp] );
    }
    
#if DEBUG_VLBI >=2
    cerr << "--- ivg::Trf Sinex::get_tropo_delays( reftype type )"
         << " : " <<tim.toc()<<" s " << endl;
#endif  
    
    return tropo;
}

// ...........................................................................
map< string, ivg::Matrix > Sinex::get_clocks( std::string & type, int & max_order )
// ...........................................................................
{
#if DEBUG_VLBI >=2
    cerr<<"+++ ivg::Trf Sinex::get_clocks( std::string )" <<endl;
    tictoc tim;
    tim.tic();
#endif
    map< string, ivg::Matrix > clo;
    
    map< string, std::vector<double> > clocks;
    map< string, std::vector<double> > clocks_std;
    map< string, std::vector<double> > epochs;
       
    // go through ALL parameter
    int ctr = 0;
    bool poly = false;
    std::vector<int> orders;
    for(vector<ivg::Param>::iterator param_iter = _parameter.begin(); param_iter!= _parameter.end(); ++param_iter)
    {         
        if( param_iter->get_type() == ivg::paramtype::clo )
        {                       
            clocks[ param_iter->get_name() ].push_back( param_iter->get_apriori()+param_iter->get_estimate() );
            clocks_std[ param_iter->get_name() ].push_back( param_iter->get_standard_deviation() );               
            epochs[ param_iter->get_name() ].push_back( param_iter->get_epoch().get_double_mjd() );
            
            orders.push_back( param_iter->get_order() );
        }       
    }
    
    // if there are not any clock parameters return empty map
    if(orders.size() == 0){
        return clo;
    }
     
    int count_offsets = std::count (orders.begin(), orders.end(), 0);
    max_order = *( std::max_element(orders.begin(), orders.end()) );
    int n = clocks.size();
    count_offsets = (double)count_offsets / (double)n;
       
    if( count_offsets == 1 )
        type = "polynomial";    
    else if( count_offsets > 1 )
        type = "cpwlf";
    else
        throw runtime_error("map< string, ivg::Matrix > Sinex::get_clocks(std::string): Unexpected parametrization of clock parameters");        
    
    
    typedef std::map<std::string, std::vector<double> >::iterator mapIter; 
    for (mapIter mit = epochs.begin(); mit != epochs.end(); mit++)
    {
        std::vector<double> vec = epochs[ mit->first ];
        vec.insert( vec.end(), clocks[ mit->first ].begin(), clocks[ mit->first ].end() );    
        vec.insert( vec.end(), clocks_std[ mit->first ].begin(), clocks_std[ mit->first ].end() );       
        ivg::Matrix T( vec.begin(), vec.end(), epochs[ mit->first ].size(), 3 );      
        clo[mit->first] = T;
    }
    
    // change station names from CDP (which is used in SINEX) to IVS names (standard in ivg::ASCOT)  
    for( auto& it : clo )
    {
        ivg::Matrix tmp = clo[ it.first ];
        clo[ _sta_assignment[it.first][ivg::staname::ivs_name] ] = tmp;
        clo.erase( _sta_assignment[it.first][ivg::staname::cdp] );
    }
    
#if DEBUG_VLBI >=2
    cerr<<"--- ivg::Trf Sinex::get_clocks( std::string )"
            <<" : "<<tim.toc()<<" s "<<endl;
#endif  
    
    return clo;
}

 std::set<std::string> Sinex::get_ref_clock_station( ){
    
    std::set<std::string> ref_stations(_station_names.begin(), _station_names.end());
    for(vector<ivg::Param>::iterator param_iter = _parameter.begin(); param_iter!= _parameter.end(); ++param_iter)
    {         
        if( param_iter->get_type() == ivg::paramtype::clo )
        {                       
            std::string station_name = _sta_assignment[ param_iter->get_name() ][ivg::staname::ivs_name];
            ref_stations.erase(station_name);
        }       
    }
    
    return ref_stations;
    
}

// ...........................................................................
int Sinex::get_num_objects(objtype type)
// ...........................................................................
{
    struct selection{vector<ivg::paramtype> types; vector<int> orders;};
    
    selection selected;
    if(type == objtype::trf)
        selected = { {ivg::paramtype::stax}, {0} };
    else if(type == objtype::eop)
        selected = { {ivg::paramtype::ut1, ivg::paramtype::xpo, ivg::paramtype::ypo, ivg::paramtype::nutx, ivg::paramtype::nuty}, {0,1} };
    else if(type == objtype::crf)
        selected = { {ivg::paramtype::ra}, {0} };
    
    int cnt = 0;
    // go through the parameter list
    for(auto &param_iter: _parameter)
    {
        // check on which parameter should be focused on, e.g. sources with ra and dec
        if(param_iter.is_type(selected.types, selected.orders))
            cnt++;
            
    }
    
    return cnt;
}
// ...........................................................................
string Sinex::_DtoE(string str)
// ...........................................................................
{
    size_t posi = str.find("D");
    if (posi!=string::npos)
    {
        str.replace(posi,1,"e");
    }

    return (str);
}
// ...........................................................................
double Sinex::_stod(string str)
// ...........................................................................
{
    double d;
    stringstream stream(str);
    stream >> d;
    return d;
}
}
