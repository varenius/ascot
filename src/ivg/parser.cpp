#include "parser.h"

namespace ivg
{
namespace parser
{

// ...........................................................................
void optl(ivg::Trf * trf_ptr,  const string path)
{
#if DEBUG_REFFRAME >=2
   cerr << "+++ void optl(ivg::Trf *,  const string )" << endl; 
   tictoc tim;
   tim.tic();
#endif
    // ...........................................................................
    ifstream inStream;
    string line;

    while( ivg::parser::get_line(path,inStream, line))
    {
        if(line.substr(0,2) != "$$")
        {
            map<string,string> tmp_names;

            string ivs_name = remove_spaces_end(line.substr( 1,8 ));
            replace_string_in_place(ivs_name, " ", "_" );

            ivg::Analysis_station * station;
            if(trf_ptr->get_station(&station, ivs_name, ivg::staname::lettercode))
            {
                ivg::Matrix optl( 3, 2 );
                istringstream optl_line( line.substr( 34, 69 ) );
                optl_line >> optl( 0 ) >> optl( 3 ) >> optl( 1 ) >> optl( 4 ) >> optl(
                              2 ) >> optl( 5 );

                station->set_optl_coeff(optl);
            }
        }
    }
    
    inStream.close();
    
#if DEBUG_REFFRAME >=2
   cerr << "--- void optl(ivg::Trf *,  const string )" << " : " << tim.toc() << " s " << endl;
#endif 
}

// ...........................................................................
void blq(ivg::Trf * trf_ptr, const string path)
// ...........................................................................
{
#if DEBUG_REFFRAME >=2
   cerr << "+++ void blq(ivg::Trf * , const string )" << endl; 
   tictoc tim;
   tim.tic();
#endif
   
    ifstream inStream;
    string line;

    while( ivg::parser::get_line(path, inStream, line) )
    {
        if(line.substr(0,2) != "$$")
        {
            string ivs_name = remove_spaces_end(line.substr( 2,8 ));
            replace_string_in_place(ivs_name, " ", "_" );

            ivg::Analysis_station * station;
            if(trf_ptr->get_station(&station, ivs_name, ivg::staname::lettercode))
            {
                ivg::parser::get_line(path,inStream, line);
                while(line.substr(0,2) == "$$")
                    ivg::parser::get_line(path, inStream, line);

                vector<float> amp_pha(66,1.0);
                for(int i = 0; i < 6 ; i++)
                {
                    int pos = i;

                    if(i>2 && i<6)
                        pos = i + 30;

                    string token;
                    stringstream tokenizer(line);
                    while(tokenizer >> token)
                    {
                        //We need to do it like this because QApplication sucks!
                        stringstream sstr; 
                        sstr << token; 
                        sstr >> amp_pha.at(pos);
                        pos = pos + 3;
                    }
                    ivg::parser::get_line(path, inStream, line);
                }
                station->set_ocean_loading_coeff(amp_pha);
            }
        }
    }
    
    inStream.close();
    
#if DEBUG_REFFRAME >=2
   cerr << "--- void blq(ivg::Trf * , const string )" << " : " << tim.toc() << " s " << endl;
#endif
}

// ...........................................................................
void harpos(ivg::Trf * trf_ptr, const string path)
// ...........................................................................
{
#if DEBUG_REFFRAME >=2
   cerr << "+++ void blq(ivg::Trf * , const string )" << endl; 
   tictoc tim;
   tim.tic();
#endif
  
    ifstream inStream;
    string line;
    ivg::Matrix phfreq(0,9);
    vector<string> harmonics;
    while( ivg::parser::get_line(path, inStream, line) )
    {
      if (line.size()>1) {
      if (line.substr(0,2) == "H ")
	{
	  
	  string nuffs=line.substr(13,46);
	  replace_string_in_place(nuffs,"D","e");
	  
	  stringstream sstrm(nuffs);
	  ivg::Matrix tmp(1,3);
	  sstrm >> tmp(0,0) >> tmp(0,1) >> tmp(0,2);
	  harmonics.push_back(remove_spaces_end(line.substr(3,8)));
	  phfreq.append_rows(tmp);
	 
	}
      if (line.substr(0,2) == "D ")
      {
	std::string ivs_name = remove_spaces_end(line.substr(13,8));
	replace_string_in_place(ivs_name, " ", "_" );

        ivg::Analysis_station * station;
	if(trf_ptr->get_station(&station, ivs_name, ivg::staname::lettercode))
            {
	      int posid;
	      ivg::Matrix sta_coef(phfreq.rows(),9,0);
	      while ((line.size()>13)&&(remove_spaces_end(line.substr(13,8))==ivs_name))
		{
		  std::string harm=remove_spaces_end(line.substr(3,8));
		  int nr;
		  for (nr=0;nr<harmonics.size();nr++)
		    if (harmonics.at(nr)==harm)
		      break;
		  stringstream sstrm(line.substr(24,55));
		  sta_coef(nr,0)=phfreq(nr,0);
		  sta_coef(nr,1)=phfreq(nr,1);
		  sta_coef(nr,2)=phfreq(nr,2);
		  
		  sstrm >> sta_coef(nr,3) >> sta_coef(nr,4) >> sta_coef(nr,5) >> sta_coef(nr,6) >> sta_coef(nr,7) >> sta_coef(nr,8);
		  posid=inStream.tellg();
		  ivg::parser::get_line(path, inStream, line);
		}
	      inStream.seekg(posid);
	      station->set_ocean_loading_coeff_harpos(sta_coef);
	    }
      }
      }
    }
    inStream.close();
    
#if DEBUG_REFFRAME >=2
   cerr << "--- void blq(ivg::Trf * , const string )" << " : " << tim.toc() << " s " << endl;
#endif
}


  
// ...........................................................................
void ecc(ivg::Trf * trf_ptr, const string path)
// ...........................................................................
{
#if DEBUG_REFFRAME >=2
   cerr << "+++ void ecc(ivg::Trf * , const string )" << endl; 
   tictoc tim;
   tim.tic();
#endif
   
    ifstream inStream;
    string line,next_line;

    while(ivg::parser::get_line(path,inStream, line))
    {
        if(line.substr(0,1) != "$")
        {
            string cdp = line.substr( 11, 4 );

            ivg::Analysis_station * station;
            if(trf_ptr->get_station(&station, cdp))
            {
                //Reading information
                ivg::Matrix new_ecc( 3, 1 );
                istringstream ecc_line( line.substr( 52, 32 ) ); // changed from 56 to 52 to read in the correct values in all cases
                ecc_line >> new_ecc( 0 ) >> new_ecc( 1 ) >> new_ecc( 2 );

                ivg::Date tmp_date(atoi(line.substr( 17,4 ).c_str()), atoi(line.substr( 22,2 ).c_str()), atoi(line.substr( 25, 2 ).c_str()));
                station->add_eccentricity(new_ecc, tmp_date, line.substr(87, 3) );
            }
        }
    }
    
    inStream.close();
    
#if DEBUG_REFFRAME >=2
   cerr << "--- void ecc(ivg::Trf * , const string )" << " : " << tim.toc() << " s " << endl;
#endif
}
  
// ...........................................................................
void gravdef(ivg::Trf * trf_ptr, const string path, ivg::Date ep)
// ...........................................................................
{
#if DEBUG_REFFRAME >=2
   cerr << "+++ void gravdef(ivg::Trf * , const string , ivg::Date  )" << endl; 
   tictoc tim;
   tim.tic();
#endif
   
    ifstream inStream;
    string line,next_line;
    ivg::parser::get_line(path,inStream, line); // First line just contain version
    while(ivg::parser::get_line(path,inStream, line))
    {
        if(line.substr(0,1) != "#")
        {
            string ivsname = remove_spaces_end(line.substr( 0,8 ));

            ivg::Analysis_station * station;
            if(trf_ptr->get_station(&station, ivsname, ivg::staname::lettercode))
            {
	  
	      int numel,i;
	        double scale2ps;
		istringstream in_line( line.substr( 9, line.length()-9) );
                in_line >> numel >> scale2ps;
		ivg::Date start(1900,1,1),end(2300,1,1);
		//Reading information
                ivg::Matrix gdel( numel, 2 );
		for (i=0;i<numel;i++)
		  {
		     ivg::parser::get_line(path,inStream, line);
		     if (line.substr(0,5) == "EPOCH")
		       {
			 
			 start=ivg::Date(atoi(line.substr(6,4).c_str()),atoi(line.substr(10,2).c_str()),atoi(line.substr(12,2).c_str()));
			 end=ivg::Date(atoi(line.substr(17,4).c_str()),atoi(line.substr(21,2).c_str()),atoi(line.substr(23,2).c_str()));
			 ivg::parser::get_line(path,inStream, line);
		       }
		     istringstream elline(line);
		     elline >> gdel(i,0) >> gdel(i,1);
		     gdel(i,1) *= scale2ps; 
		  }

		if((ep>=start)&&(ep<=end))
		  {
		    station->set_gravdef(gdel);
		  }
            }
        }
    }
    
    inStream.close();
    
#if DEBUG_REFFRAME >=2
   cerr << "--- void gravdef(ivg::Trf * , const string , ivg::Date )" << " : " << tim.toc() << " s " << endl;
#endif
}
  
// ...........................................................................
void bindisp(ivg::Trf * trf_ptr, const string folder_path, map<string,string> correspondence, ivg::Date start, ivg::Date end)
// ...........................................................................
{
#if DEBUG_REFFRAME >=2
   cerr << "+++ void bindisp(ivg::Trf * , const string , ivg::Date , ivg::Date )" << endl; 
   tictoc tim;
   tim.tic();
#endif
   
    double start_mjd = start.get_double_mjd();
    double end_mjd = end.get_double_mjd();

    vector<string> station_names = trf_ptr->get_station_names(ivg::staname::ivs_name);

    for (int i=0 ; i<station_names.size(); i++)
    {
        ivg::Analysis_station * station;
        trf_ptr->get_station(&station, station_names.at(i), ivg::staname::lettercode);

        ifstream infile;
        string path = folder_path + station->get_name(ivg::staname::ivs_name) +".bds";

        infile.open(path, ios::binary | ios::in);

        if( !infile.is_open() )
        {
//            log<DETAIL>("*** void bindisp: ") % path % " not existent. " % station->get_name(ivg::staname::corres) % ".bds will be used instead.";
            path = folder_path + correspondence[station->get_name(ivg::staname::ivs_name)] + ".bds";
            infile.open(path, ios::binary | ios::in);
        }

        // if the file doesn't exist, we skip the analysis center
        if(!file_exists(path))
        {
            log<WARNING>("!!! Parsing NTAPL data: ") % path % " does not exist. Skipping " %  station->get_name(ivg::staname::ivs_name);
            continue;
        }
        
        char tmp;
        stringstream ss;
        for(int i=0; i<8 ; i++)
        {
            infile.read((char*)&tmp, sizeof(tmp));
            ss << tmp;
        }

        if(ss.str() == "BINDISP ")
        {
            int revision_mjd;
            infile.read((char*)&revision_mjd, sizeof(revision_mjd));

            char int_format,float_format;
            infile.read((char*)&int_format, sizeof(int_format));
            infile.read((char*)&float_format, sizeof(float_format));

            short reserved;
            infile.read((char*)&reserved, sizeof(reserved));

            ss.str(string());
            for(int i=0; i<8 ; i++)
            {
                infile.read((char*)&tmp, sizeof(tmp));
                ss << tmp;
            }
            string site = ss.str();

            int n_records;
            infile.read((char*)&n_records, sizeof(n_records));

            float interval;
            infile.read((char*)&interval, sizeof(interval));

            ivg::Matrix xyz0( 3, 1 );
            infile.read((char*)&xyz0(0,0), sizeof(double));
            infile.read((char*)&xyz0(1,0), sizeof(double));
            infile.read((char*)&xyz0(2,0), sizeof(double));

            int initial_mjd;
            infile.read((char*)&initial_mjd, sizeof(initial_mjd));
            float initial_secs;
            infile.read((char*)&initial_secs, sizeof(initial_secs));

//                    cout << " | Info-Rows BINDISP " << site << " |" << endl;
//                    cout << " | revision_mjd: " << revision_mjd << " | int_format: " << int_format << " | float_format: " << float_format;
//                    cout << " | n_records: " << n_records << " | interval: " << interval << " |" << endl;
//                    cout << " | xyz0: " << xyz0(0,0) << " " << xyz0(1,0) << " " << xyz0(2,0);
//                    cout << " | initial_mjd: " << initial_mjd << " | initial_secs: " << initial_secs << " |" << endl;

            int elements=0;
            short x,y,z;
            vector<double> exyz_vec; //Epoch, x, y, z
            double interval_mjd = ((double)interval / 86400.0); //0.25 Tage
            double epoch_mjd = (double)initial_mjd + ((double)initial_secs / 86400.0);

            for(int j=0 ; j<n_records ; j++)
            {
                infile.read((char*)&x, sizeof(x));
                infile.read((char*)&y, sizeof(y));
                infile.read((char*)&z, sizeof(z));
                infile.read((char*)&reserved, sizeof(reserved));

                if(epoch_mjd >= start_mjd && epoch_mjd <= end_mjd)
                {
                    exyz_vec.push_back(epoch_mjd);
                    exyz_vec.push_back((double) x * 1.0e-5);
                    exyz_vec.push_back((double) y * 1.0e-5);
                    exyz_vec.push_back((double) z * 1.0e-5);
                    elements++;
                }
                epoch_mjd = epoch_mjd + interval_mjd;
            }

            ivg::Matrix exyz(exyz_vec);
            exyz.resize(4,elements);

            station->set_nontidal_aplo(exyz);
        }
        infile.close();
    }
#if DEBUG_REFFRAME >=2
   cerr << "--- void bindisp(ivg::Trf * , const string , ivg::Date , ivg::Date )" << " : " << tim.toc() << " s " << endl;
#endif
}

  // ...........................................................................
void vie_atm(ivg::Trf * trf_ptr, const string folder_path, map<string,string> correspondence, ivg::Date start, ivg::Date end)
// ...........................................................................
{
#if DEBUG_REFFRAME >=2
   cerr << "+++ void bindisp(ivg::Trf * , const string , ivg::Date , ivg::Date )" << endl; 
   tictoc tim;
   tim.tic();
#endif
   
    double start_mjd = start.get_double_mjd();
    double end_mjd = end.get_double_mjd();

    vector<string> station_names = trf_ptr->get_station_names(ivg::staname::ivs_name);
    int yr1 = start.get_int_year();
    int yr2 = end.get_int_year();
    vector<vector<double>> exyz_vec; //Epoch, x, y, z
    vector<double> du,de,dn;
    for (int i=0 ; i<station_names.size(); i++)
      {

	 ivg::Analysis_station * station;
         trf_ptr->get_station(&station, station_names.at(i), ivg::staname::lettercode);
	 ivg::Matrix llh=station->calc_lat_lon_h();
	 du.push_back(cos(llh(1))*cos(llh(0)));
	 du.push_back(sin(llh(1))*cos(llh(0)));
	 du.push_back(sin(llh(0)));
	 dn.push_back(-cos(llh(1))*sin(llh(0)));
	 dn.push_back(-sin(llh(1))*sin(llh(0)));
	 dn.push_back(cos(llh(0)));
	 de.push_back(-sin(llh(1)));
	 de.push_back(cos(llh(1)));
	 de.push_back(0.0);
	 exyz_vec.push_back({});
	


      }
    
    for (int yr=yr1;yr<=yr2;yr++)
      {
	
	string path = folder_path + "/y" + std::to_string(yr) + ".apl_r";
	
	if (file_exists(path))
	  {
	    
	    ifstream infile(path);
	    string line,sta;
	    double mjd,up,east,north;
	    while (getline(infile,line))
	      {
		if (line.substr(0,1)!="!")
		  {
		    stringstream ss(line);
		    ss >> sta >> mjd >> up >> east >> north;
		    if ((mjd>=start_mjd)&&(mjd<=end_mjd)) {
		      for (int i=0 ; i<station_names.size(); i++) {
			
			if ((station_names[i]==sta)||(correspondence[station_names[i]]==sta)){
			 
			  exyz_vec[i].push_back(mjd);
			  
			  exyz_vec[i].push_back(up*du[i*3]+east*de[i*3]+north*dn[i*3]);
			  
			  exyz_vec[i].push_back(up*du[i*3+1]+east*de[i*3+1]+north*dn[i*3+1]);
			  
			  exyz_vec[i].push_back(up*du[i*3+2]+east*de[i*3+2]+north*dn[i*3+2]);
			 
			  
			}
			  
		      }

		    }
		  }
		
	      }
	    infile.close();
	  }
	else
	  cerr << "Atmp loading file " << path << " does not exist" <<endl;
      }
    for (int i=0 ; i<station_names.size(); i++)
      {
	 ivg::Analysis_station * station;
         trf_ptr->get_station(&station, station_names.at(i), ivg::staname::lettercode);
	 ivg::Matrix exyz(exyz_vec[i]);
         exyz.resize(4,int(exyz_vec[i].size()/4));
	 for (int j=1;j<exyz.cols();j++) {
	   if (exyz(0,j)<=exyz(0,j-1))
	     exyz.rem_c(j);
	 }
	 station->set_nontidal_aplo(exyz);
      }
    
  
#if DEBUG_REFFRAME >=2
   cerr << "--- void vie_atm(ivg::Trf * , const string , ivg::Date , ivg::Date )" << " : " << tim.toc() << " s " << endl;
#endif
}

// ...........................................................................
void hps(ivg::Trf * trf_ptr, const string path)
// ...........................................................................
{
#if DEBUG_REFFRAME >=2
   cerr << "+++ void hps(ivg::Trf * , const string )" << endl; 
   tictoc tim;
   tim.tic();
#endif
   
    map<string,ivg::Wave> waves;
    
    map<string,string> assignment;

    map<string,ivg::Matrix> XYZ;
    
    ifstream inStream;
    string line;

    bool get_inside=true;
    while( ivg::parser::get_line(path, inStream, line) )
    {
        
        if(line.substr(0,1) == "S")
        {
            ivg::Matrix xyz( 3, 1 );
            istringstream xyz_line( line.substr( 13, 41 ) );
            xyz_line >> xyz( 0 ) >> xyz( 1 ) >> xyz( 2 );
            
            XYZ[remove_spaces_end(line.substr( 3, 8 ))] = xyz;
            
            
        }
        else if(line.substr(0,1) == "H" && line.substr(0,6) != "HARPOS")
        {
            ivg::Wave tmp;
            tmp.name = remove_spaces_end(line.substr( 3, 8 ));
            string str = line.substr( 13, 13 );
            replace_string_in_place( str, "D", "e" );
            tmp.phase = s2d( str );

            str = line.substr( 28,19 );
            replace_string_in_place( str, "D", "e" );
            tmp.freq = s2d( str );

            str = line.substr( 49,10 );
            replace_string_in_place( str, "D", "e" );
            tmp.accel = s2d( str );

            waves[tmp.name] = tmp;

        }
        else if(line.substr(0,1) == "D")
        {
            //Only get inside the first time, to detect the correct assignment
            if(get_inside)
            {    
                for (std::vector<ivg::Analysis_station>::iterator sta_iter = trf_ptr->begin(); sta_iter!= trf_ptr->end(); sta_iter++)
                {
                    double norm = 999999999;
                    for(map<string, ivg::Matrix>::iterator data_iter = XYZ.begin(); data_iter != XYZ.end(); ++data_iter )
                    {
                       ivg::Matrix norm_matrix = (sta_iter->get_xyz0().get_sub(0,0,2,0) - data_iter->second).norm();
                       
                       if(norm_matrix(0) < norm){
                           norm = norm_matrix(0);
                           
                           //if distance less than 3000m between sites, take the data
                           if(norm <= 3000.0)
                               assignment[sta_iter->get_name(ivg::staname::ivs_name)] = data_iter->first;
                       }
                    }
                }
            get_inside=false;
            }
            
            string wname = remove_spaces_end(line.substr( 3, 8 ));
            string corres_name = remove_spaces_end(line.substr( 13, 8 ));

            for (std::vector<ivg::Analysis_station>::iterator sta_iter = trf_ptr->begin(); sta_iter!= trf_ptr->end(); sta_iter++)
            {
                if(assignment[sta_iter->get_name(ivg::staname::ivs_name)] == corres_name)
                {
                    istringstream cos_sin_line( line.substr( 24, 55 ) );
                    cos_sin_line >> waves[wname].up_cos >> waves[wname].east_cos >>
                                 waves[wname].north_cos;
                    cos_sin_line >> waves[wname].up_sin >> waves[wname].east_sin >>
                                 waves[wname].north_sin;

                sta_iter->set_tidal_aplo(waves);
                }
            }
        }
    }
    
    inStream.close();
    
#if DEBUG_REFFRAME >=2
   cerr << "--- void hps(ivg::Trf * , const string )" << " : " << tim.toc() << " s " << endl;
#endif
}
// ...........................................................................
void dat(ivg::Trf * trf_ptr, const string path)
// ...........................................................................
{
#if DEBUG_REFFRAME >=2
   cerr << "+++ void dat(ivg::Trf * , const string )" << endl; 
   tictoc tim;
   tim.tic();
#endif
    // THIS FUNCTION IS NOT REALLY VALIDATED
    // IDEA: For Tidal Athmospheric Pressure Loading the hps file does not contain all requested stations.
    // Using *.dat files from http://ggosatm.hg.tuwien.ac.at/LOADING/VERSION4/TIDAL/ could solve this issue.
    // Effect of TAPL very small and not within conventions, so not such a high priority.
   
    map<string,ivg::Wave> waves;
    
    ifstream inStream;
    string line;

    bool get_inside=true;
    while( ivg::parser::get_line(path, inStream, line) )
    {
        string ivs_name = remove_spaces_end(line.substr(0,8));
        
        ivg::Analysis_station * station;
        if(trf_ptr->get_station(&station, ivs_name, ivg::staname::lettercode))
        {
            vector<string> tokens = get_tokens(line.substr(11));

            ivg::Wave S1;
            S1.name = "S1";
            S1.accel = 0.0;
            S1.phase = M_PI;
            S1.freq = 7.272205216643e-05;

            ivg::Wave S2;
            S2.name = "S2";
            S2.accel = 0.0;
            S2.phase = 0.0;
            S2.freq = 1.454441043329e-04;

            S1.up_cos = s2d(tokens.at(0))/1000.0;
            S1.up_sin = s2d(tokens.at(1))/1000.0;
            S2.up_cos = s2d(tokens.at(2))/1000.0;
            S2.up_sin = s2d(tokens.at(3))/1000.0;

            S1.east_cos = s2d(tokens.at(6))/1000.0;
            S1.east_sin = s2d(tokens.at(7))/1000.0;
            S2.east_cos = s2d(tokens.at(8))/1000.0;
            S2.east_sin = s2d(tokens.at(9))/1000.0;

            S1.north_cos = s2d(tokens.at(12))/1000.0;
            S1.north_sin = s2d(tokens.at(13))/1000.0;
            S2.north_cos = s2d(tokens.at(14))/1000.0;
            S2.north_sin = s2d(tokens.at(15))/1000.0;
	   
            waves[S1.name] = S1;
            waves[S2.name] = S2;
            
            station->set_tidal_aplo(waves);
        }
    }
    
    inStream.close();
    
#if DEBUG_REFFRAME >=2
   cerr << "--- void dat(ivg::Trf * , const string )" << " : " << tim.toc() << " s " << endl;
#endif
}
// ...........................................................................
void antenna_info(ivg::Trf * trf_ptr,  const string path)
{
// ...........................................................................
#if DEBUG_REFFRAME >=2
   cerr << "+++ void antenna_info(ivg::Trf * ,  const string )" << endl; 
   tictoc tim;
   tim.tic();
#endif
   
    ifstream inStream;
    string line;

    while( ivg::parser::get_line(path,inStream, line))
    {
        if(line.substr(0,12) == "ANTENNA_INFO")
        {
            string ivs_name = remove_spaces_end(line.substr( 14, 8 ));

            ivg::Analysis_station * station;
            if(trf_ptr->get_station(&station, ivs_name, ivg::staname::lettercode))
            {
                ivg::Antenna tmp =
                {
                    remove_spaces_end(line.substr( 24, 7 )),remove_spaces_end(line.substr( 32, 7 )),
                    remove_spaces_end(line.substr( 40, 7 )),remove_spaces_end(line.substr( 48, 7 )),
                    s2d(line.substr( 57, 4 )),s2d(line.substr( 62, 4 )),s2d(line.substr( 67, 4 )),
                    s2d(line.substr( 72, 6 )),s2d(line.substr( 80, 5 )),s2d(line.substr( 86, 7 )),
                    s2d(line.substr( 94, 6 )),s2d(line.substr( 102, 7 )),s2d(line.substr( 111, 7 )),
                    s2d(line.substr( 119, 7 )),s2d(line.substr( 128, 7 )),s2d(line.substr( 136, 7 )),
                    s2d(line.substr( 145, 7 )),s2d(line.substr( 153, 7 )),s2d(line.substr( 162, 7 )),
                    s2d(line.substr( 170, 7 )),
                };

                station->set_antenna_info(tmp);
            }
        }
    }
    
    inStream.close();
    
#if DEBUG_REFFRAME >=2
   cerr << "--- void antenna_info(ivg::Trf * ,  const string )" << " : " << tim.toc() << " s " << endl;
#endif
}

// ...........................................................................
void external_met_data( ivg::Trf * trf_ptr, const string path, ivg::Date start,
                        ivg::Date end, ivg::extdata type )
// ...........................................................................
{
#if DEBUG_REFFRAME >=2
   cerr << "+++ void external_met_data(ivg::Trf * , const string , ivg::Date , ivg::Date )" << endl; 
   tictoc tim;
   tim.tic();
#endif
   
    double start_mjd = start.get_double_mjd();
    double end_mjd = end.get_double_mjd();

    vector<string> station_names = trf_ptr->get_station_names(ivg::staname::ivs_name);

    for (int i=0 ; i<station_names.size(); i++)
    {
        ivg::Analysis_station * station;
        if(trf_ptr->get_station(&station, station_names.at(i), ivg::staname::lettercode))
        {
            string load_path = path + station->get_name(ivg::staname::ivs_name) + ".bin";
            
            // if the file doesnt exist, we skip the analysis center
            if(!file_exists(load_path))
            {
                log<WARNING>("!!! Parsing VMF1 data: ") % load_path % " does not exist. Skipping " %  station->get_name(ivg::staname::ivs_name);
                continue;
            }
            
            ivg::Matrix vmf;
            vmf.load_bin(load_path);

            int start_pos  = (vmf.get_col(0) - start_mjd).abs().minIdx();
            int end_pos  = (vmf.get_col(0) - end_mjd).abs().minIdx();
//                  cout << "Start MJD: " << start_mjd << " / Ende MJD: " << end_mjd << " / mini_index: " << start_pos << " / maxi_index: " << end_pos << " / gesamtlang: " << vmf1.rows() << endl;
            
            if( type == ivg::extdata::MAPPING )
                station->set_vmf1_data(vmf.get_sub(start_pos,0,end_pos,9));
	    else if( type == ivg::extdata::MAPPING3 )
                station->set_vmf3_data(vmf.get_sub(start_pos,0,end_pos,9));
            else if( type == ivg::extdata::HYDROSTATIC )
                station->set_zhd_data(vmf.get_sub(start_pos,0,end_pos,9));
        }
    }
#if DEBUG_REFFRAME >=2
   cerr << "--- void external_met_data(ivg::Trf * , const string , ivg::Date , ivg::Date )" << " : " << tim.toc() << " s " << endl;
#endif
}
// ...........................................................................
void external_met_data_ascii( ivg::Trf * trf_ptr, const string path, ivg::Date start,
                        ivg::Date end, ivg::extdata type )
// ...........................................................................
{
#if DEBUG_REFFRAME >=2
   cerr << "+++ void external_met_data_ascii(ivg::Trf * , const string , ivg::Date , ivg::Date )" << endl; 
   tictoc tim;
   tim.tic();
#endif
   
    double start_mjd = start.get_double_mjd();
    double end_mjd = end.get_double_mjd();

    vector<string> station_names = trf_ptr->get_station_names(ivg::staname::ivs_name);
    vector<ivg::Matrix> vmf(station_names.size(),ivg::Matrix(0,10));
    
    int yr1 = start.get_int_year();
    int yr2 = end.get_int_year();
    
    for (int yr=yr1;yr<=yr2;yr++)
      {
	string vmffile;
	if (type ==  ivg::extdata::MAPPING )
	  vmffile=path + "/y" + std::to_string(yr) + ".vmf1_r";
	else if (type ==  ivg::extdata::MAPPING3 )
	  vmffile=path + "/y" + std::to_string(yr) + ".vmf3_r";
	else
	  vmffile=path + "/y" + std::to_string(yr) + ".vmf1_r";

        if (file_exists(vmffile))
          {

            ifstream infile(vmffile);
            string line,sta;
	    double mjd;
	     while (getline(infile,line))
              {
                if (line.substr(0,1)!="#")
		  {
		    stringstream ss(line);
		    ss >> sta;
		    ss >> mjd;
		    if ((mjd>=start_mjd)&&(mjd<=end_mjd)) {
                      for (int i=0 ; i<station_names.size(); i++) {
			
                        if (station_names[i]==sta)
			  {
			   
			    ivg::Matrix row(1,10);
			    row(0,0)=mjd;
			    for (int j=1;j<6;j++)
			      ss >> row(0,j);
			    row(0,6)=0;
			    for (int j=7;j<8;j++)
			      ss >> row(0,j);
			    row(0,9)=0;
			    vmf[i].append_rows(row);
			  }
		      }
		    }
		  }
	      }
	  }
	else
	  cerr << "File does not exist: " << vmffile << endl;
      }
    
    
    for (int i=0 ; i<station_names.size(); i++)
    {
        ivg::Analysis_station * station;
        if(trf_ptr->get_station(&station, station_names.at(i), ivg::staname::lettercode))
        {
            
            
            if( type == ivg::extdata::MAPPING )
                station->set_vmf1_data(vmf[i]);
	    else if( type == ivg::extdata::MAPPING3 )
                station->set_vmf3_data(vmf[i]);
            else if( type == ivg::extdata::HYDROSTATIC )
                station->set_zhd_data(vmf[i]);
        }
    }
#if DEBUG_REFFRAME >=2
   cerr << "--- void external_met_data_ascii(ivg::Trf * , const string , ivg::Date , ivg::Date )" << " : " << tim.toc() << " s " << endl;
#endif
}

// ...........................................................................
void hydlo(ivg::Trf * trf_ptr, const string folder_path)
// ...........................................................................
{
#if DEBUG_REFFRAME >=2
   cerr << "+++ void hydlo(ivg::Trf * , const string )" << endl; 
   tictoc tim;
   tim.tic();
#endif
   
    vector<string> station_names = trf_ptr->get_station_names(
                                       ivg::staname::ivs_name);

    for (int i=0 ; i<station_names.size(); i++)
    {
        ivg::Analysis_station * station;
        trf_ptr->get_station(&station, station_names.at(i));

        string ending;
        if (folder_path.find("cmte") != std::string::npos)
            ending = "_CMTE.txt";
        else if (folder_path.find("cmse") != std::string::npos)
            ending = "_CMSE.txt";

        string path = folder_path + station->get_name(ivg::staname::ivs_name) + ending;

        if(file_exists(path))
        {
            ifstream inStream;
            string line;

            ivg::Matrix hydlo( 1, 4 , 0.0);

            while(ivg::parser::get_line(path,inStream, line))
            {
                if(line.substr(0,1) != "#")
                {
                    ivg::Matrix mjdxyz( 1, 4 );
                    istringstream mjdxyz_line( line.substr( 0, 48 ) );
                    mjdxyz_line >> mjdxyz( 0 ) >> mjdxyz( 1 ) >> mjdxyz( 2 ) >> mjdxyz( 3 );
                    hydlo.append_rows(mjdxyz);
                }
            }
            hydlo.rem_r(0);
            station->set_hydlo(hydlo);
            inStream.close();
        }
        else
            log<WARNING>("!!! Hydrological data not available. Failed to open ") % path;
    }
#if DEBUG_REFFRAME >=2
   cerr << "--- void hydlo(ivg::Trf * , const string )" << " : " << tim.toc() << " s " << endl;
#endif
}

// ...........................................................................
void seasonals(ivg::Trf * trf_ptr, const string path)
// ...........................................................................
{
  ifstream inStream;
  string line;
  vector<double> period;
  int freqnr;
  map<string, map<int, map<string, vector<double> > > >  savemap;
  while(ivg::parser::get_line(path,inStream, line))
    {
      if(line.find("Frequency")!=string::npos)
	{
	  double tmp;
	 
	  stringstream sstr1,sstr2; 
	  sstr1 << remove_spaces_end(line.substr(10,2)); 
	  sstr1 >> freqnr;
	  sstr2 << remove_spaces_end(line.substr(15,7)); 
	  sstr2 >> tmp;
	  period.push_back(tmp);
	 
	}
      if (line.substr(0,1)==" " && (line.substr(22,1)=="X" || line.substr(22,1)=="Y" || line.substr(22,1)=="Z"))
	{
	  string domes;
	  domes=line.substr(9,9);
	  vector<double> coeff(4,0.0);
	  if (line.size()>36) {
	    stringstream tokenizer(line.substr(23,32));
            
	    double val;
	    int j=0;
	    
	    while(tokenizer >> val)
	      {
		coeff.at(j) = val;
		j++;
	      }
	  }
	  
	  savemap[domes][period.size()-1][line.substr(22,1)]=coeff;
	    
	}
    }
  inStream.close();
  for(auto &station: savemap)
    {
      ivg::Analysis_station * station_ptr;
      if(trf_ptr->get_station(&station_ptr, station.first))
	{
	  ivg::Matrix coefs;
	  for(auto &curper: station.second)
	    {
	      vector<double> tmp;
	      tmp.push_back(period[curper.first]);
	     
	      tmp.insert(tmp.end(),curper.second["X"].begin(),curper.second["X"].end());
	      
	      tmp.insert(tmp.end(),curper.second["Y"].begin(),curper.second["Y"].end());
	     
	      tmp.insert(tmp.end(),curper.second["Z"].begin(),curper.second["Z"].end());
	      
	      coefs.append_rows(ivg::Matrix(tmp).transpose());
	    
	  
	    }
	  station_ptr->set_seasonals(coefs);
	  
	}
    }
}

// ...........................................................................
void psd_coefficients(ivg::Trf * trf_ptr, const string path)
// ...........................................................................
{
#if DEBUG_REFFRAME >=2
   cerr << "+++ void psd_coefficients(ivg::Trf * trf_ptr, const string file_path)" << endl; 
   tictoc tim;
   tim.tic();
#endif
   
    string line;
    ifstream inStream;
    
    map< ivg::Analysis_station*, vector<Psd> > ac_psd_assignment;
    
     // in case of snx-psd-file
    if(path.find("snx") != string::npos)
    {
        map<string, map<double, map<string, map< string, vector<double> > > > > savemap;
        
        while(ivg::parser::get_line(path,inStream, line))
        {
            if(line.find("+SOLUTION/ESTIMATE")!=string::npos)
            {
                while (getline(inStream,line,'\n') && line.substr(0,1)!="-")
                {
                    if(line.substr(0,1) != "*")
                    {
                        string type = remove_spaces_end(line.substr(7,6));
                        string code = line.substr(14,4).c_str();
                        ivg::Date epoch(line.substr(27,12),"SINEX");
                        double mjd = epoch.get_double_mjd();
                        string component = type.substr(5,1);
                        string mode = type.substr(1,3);
                        
                        double estimate;

                        stringstream sstr1; 
                        sstr1 << remove_spaces_end(line.substr(47,21)); 
                        sstr1 >> estimate;
                        
                        if(type.substr(0,1) == "A")
                            estimate *= 1e3;

                        savemap[code][mjd][component][mode].push_back(estimate);
                    }
                } 
            }
        }
        
        
        for(auto &station: savemap)
        {
            for(auto &mjd: station.second)
            {
//                cerr << "station: " << station.first << " / mjd: " << mjd.first << endl;
                
                ivg::Psd tmp_psd = { mjd.first, 0, 0, 0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
                
                for(auto &comp: mjd.second)
                {
                    // E, H, N
                    for(auto &mode: comp.second)
                    {
                        // EXP, LOG
                        int n_modes = mode.second.size();
                        if(mode.first == "EXP" && comp.first == "E")
                            tmp_psd.e_mode += n_modes;
                        else if(mode.first == "LOG" && comp.first == "E")
                            tmp_psd.e_mode += 2*n_modes -3;
                        else if(mode.first == "EXP" && comp.first == "N" )
                            tmp_psd.n_mode += n_modes;
                        else if(mode.first == "LOG" && comp.first == "N")
                            tmp_psd.n_mode += 2*n_modes -3;
                        else if(mode.first == "EXP" && comp.first == "H" )
                            tmp_psd.u_mode += n_modes;
                        else if(mode.first == "LOG" && comp.first == "H")
                            tmp_psd.u_mode += 2*n_modes -3;
                    }                    
                }
                
                // we got five different modes
                if(tmp_psd.e_mode == 1 )
                {
                    tmp_psd.e_a1 = mjd.second["E"]["LOG"].at(0);
                    tmp_psd.e_t1 = mjd.second["E"]["LOG"].at(1);
                }
                else if(tmp_psd.e_mode == 2)
                {
                    tmp_psd.e_a1 = mjd.second["E"]["EXP"].at(0);
                    tmp_psd.e_t1 = mjd.second["E"]["EXP"].at(1);
                }
                else if(tmp_psd.e_mode == 3)
                {
                    tmp_psd.e_a1 = mjd.second["E"]["LOG"].at(0);
                    tmp_psd.e_t1 = mjd.second["E"]["LOG"].at(1);
                    tmp_psd.e_a2 = mjd.second["E"]["EXP"].at(0);
                    tmp_psd.e_t2 = mjd.second["E"]["EXP"].at(1);
                }
                else if(tmp_psd.e_mode == 4)
                {
                    tmp_psd.e_a1 = mjd.second["E"]["EXP"].at(0);
                    tmp_psd.e_t1 = mjd.second["E"]["EXP"].at(1);
                    tmp_psd.e_a2 = mjd.second["E"]["EXP"].at(2);
                    tmp_psd.e_t2 = mjd.second["E"]["EXP"].at(3);
                }
		else if(tmp_psd.e_mode == 5)
                {
                    tmp_psd.e_a1 = mjd.second["E"]["LOG"].at(0);
                    tmp_psd.e_t1 = mjd.second["E"]["LOG"].at(1);
                    tmp_psd.e_a2 = mjd.second["E"]["LOG"].at(2);
                    tmp_psd.e_t2 = mjd.second["E"]["LOG"].at(3);
                }
                
                // we got four different modes
                if(tmp_psd.n_mode == 1 )
                {
                    tmp_psd.n_a1 = mjd.second["N"]["LOG"].at(0);
                    tmp_psd.n_t1 = mjd.second["N"]["LOG"].at(1);
                }
                else if(tmp_psd.n_mode == 2)
                {
                    tmp_psd.n_a1 = mjd.second["N"]["EXP"].at(0);
                    tmp_psd.n_t1 = mjd.second["N"]["EXP"].at(1);
                }
                else if(tmp_psd.n_mode == 3)
                {
                    tmp_psd.n_a1 = mjd.second["N"]["LOG"].at(0);
                    tmp_psd.n_t1 = mjd.second["N"]["LOG"].at(1);
                    tmp_psd.n_a2 = mjd.second["N"]["EXP"].at(0);
                    tmp_psd.n_t2 = mjd.second["N"]["EXP"].at(1);
                }
                else if(tmp_psd.n_mode == 4)
                {
                    tmp_psd.n_a1 = mjd.second["N"]["EXP"].at(0);
                    tmp_psd.n_t1 = mjd.second["N"]["EXP"].at(1);
                    tmp_psd.n_a2 = mjd.second["N"]["EXP"].at(2);
                    tmp_psd.n_t2 = mjd.second["N"]["EXP"].at(3);
                }
                else if(tmp_psd.n_mode == 5)
                {
                    tmp_psd.n_a1 = mjd.second["N"]["LOG"].at(0);
                    tmp_psd.n_t1 = mjd.second["N"]["LOG"].at(1);
                    tmp_psd.n_a2 = mjd.second["N"]["LOG"].at(2);
                    tmp_psd.n_t2 = mjd.second["N"]["LOG"].at(3);
                }
                // we got five different modes
                if(tmp_psd.u_mode == 1 )
                {
                    tmp_psd.u_a1 = mjd.second["H"]["LOG"].at(0);
                    tmp_psd.u_t1 = mjd.second["H"]["LOG"].at(1);
                }
                else if(tmp_psd.u_mode == 2)
                {
                    tmp_psd.u_a1 = mjd.second["H"]["EXP"].at(0);
                    tmp_psd.u_t1 = mjd.second["H"]["EXP"].at(1);
                }
                else if(tmp_psd.u_mode == 3)
                {
                    tmp_psd.u_a1 = mjd.second["H"]["LOG"].at(0);
                    tmp_psd.u_t1 = mjd.second["H"]["LOG"].at(1);
                    tmp_psd.u_a2 = mjd.second["H"]["EXP"].at(0);
                    tmp_psd.u_t2 = mjd.second["H"]["EXP"].at(1);
                }
                else if(tmp_psd.u_mode == 4)
                {
                    tmp_psd.u_a1 = mjd.second["H"]["EXP"].at(0);
                    tmp_psd.u_t1 = mjd.second["H"]["EXP"].at(1);
                    tmp_psd.u_a2 = mjd.second["H"]["EXP"].at(2);
                    tmp_psd.u_t2 = mjd.second["H"]["EXP"].at(3);
                }
                else if(tmp_psd.u_mode == 5)
                {
                    tmp_psd.u_a1 = mjd.second["H"]["LOG"].at(0);
                    tmp_psd.u_t1 = mjd.second["H"]["LOG"].at(1);
                    tmp_psd.u_a2 = mjd.second["H"]["LOG"].at(2);
                    tmp_psd.u_t2 = mjd.second["H"]["LOG"].at(3);
                }
                
//                cerr << "E_mode " << tmp_psd.e_mode << ": " << tmp_psd.e_a1 << "|" << tmp_psd.e_t1 << "|" << tmp_psd.e_a2 << "|" << tmp_psd.e_t2 << endl;
//                cerr << "N_mode " << tmp_psd.n_mode << ": " << tmp_psd.n_a1 << "|" << tmp_psd.n_t1 << "|" << tmp_psd.n_a2 << "|" << tmp_psd.n_t2 << endl;
//                cerr << "U_mode " << tmp_psd.u_mode << ": " << tmp_psd.u_a1 << "|" << tmp_psd.u_t1 << "|" << tmp_psd.u_a2 << "|" << tmp_psd.u_t2 << endl;
                
                ivg::Analysis_station * station_ptr;
                if(trf_ptr->get_station(&station_ptr, station.first))
                {
                    ac_psd_assignment[station_ptr].push_back(tmp_psd);
                    ivg::Date tmp_date(tmp_psd.mjd);
                    log<INFO>("*** Loading psd-coefficients for station ") % station_ptr->get_name(ivg::staname::ivs_name) % " at epoch: " % tmp_date.get_date_time("YY:DOY:SSSSS");
                } 
            }
        }
        
    }
    else
    {

        while(ivg::parser::get_line(path,inStream, line))
        {
            // get information about the site
            string cdp   = line.substr(1,4);
            string domes = line.substr(9,9);

            ivg::Analysis_station * station;
            if(trf_ptr->get_station(&station, domes))
            {
                //vector<Psd> tmp_vec;
                //ac_psd_assignment[station] = tmp_vec;

                // post-seismic deformations (coefficients)
                Psd tmp;
                tmp.mjd = ivg::Date(line.substr(19,12), "SINEX" ).get_double_mjd();

                log<INFO>("*** Loading psd-coefficients for station ") % station->get_name(ivg::staname::ivs_name) % " at epoch: " % line.substr(19,12);

                double dpos;
                ivg::Matrix d_latlonh(3,1,0.0);
                for(int i=0;i<3;++i)
                {
                    if(i>0)
                        ivg::parser::get_line(path,inStream, line);
		   
		    vector<double> coeff(4,0.0);
		    if (line.size()>36) {
		      stringstream tokenizer(line.substr(36,31));
                   
		      double val;
		      int j=0;
		      while(tokenizer >> val)
                        coeff.at(j++) = val;
		    }


                    string comp = line.substr(32,1);
                    int mode = stoi(line.substr(34,1));

                    if(comp == "E")
                    {
                        tmp.e_mode = mode;
                        tmp.e_a1 = coeff.at(0);
                        tmp.e_t1 = coeff.at(1);
                        tmp.e_a2 = coeff.at(2);
                        tmp.e_t2 = coeff.at(3);
                    }
                    else if( comp == "N" )
                    {
                        tmp.n_mode = mode;
                        tmp.n_a1 = coeff.at(0);
                        tmp.n_t1 = coeff.at(1);
                        tmp.n_a2 = coeff.at(2);
                        tmp.n_t2 = coeff.at(3);
                    }
                    else if( comp == "U" )
                    {
                        tmp.u_mode = mode;
                        tmp.u_a1 = coeff.at(0);
                        tmp.u_t1 = coeff.at(1);
                        tmp.u_a2 = coeff.at(2);
                        tmp.u_t2 = coeff.at(3);
                    }
                }
                ac_psd_assignment[station].push_back(tmp);
            }
        }
        
    }
    
    inStream.close();
    
    for(auto &ac: ac_psd_assignment)
        ac.first->set_psd_coeff(ac.second);
   
   
#if DEBUG_REFFRAME >=2
    cerr << "--- void psd_coefficients(ivg::Trf * trf_ptr, const string file_path)" << " : " << tim.toc() << " s " << endl;
#endif  
}

// ...........................................................................
void raytraced_delays( ivg::Trf * trf_ptr, const string path )
{
// ...........................................................................
#if DEBUG_REFFRAME >=2
   cerr << "+++ void raytraced_delays(ivg::Trf * , const string )" << endl; 
   tictoc tim;
   tim.tic();
#endif
   
    ifstream inStream;
    string line;
    std::map< std::string,std::map< ivg::Date,std::map<std::string,ivg::Matrix> > > rtmap;
    
    while( ivg::parser::get_line(path, inStream, line))
    {       
        // read only '0-records section' 
        // (based on RADIATE raytracing program by A. Hofmeister, TU Vienna)
        if(line.substr(0,1) == "O")
        {
            // read zenith hydrostatic and wet delay as well as wet mapping factor
            ivg::Matrix atm( 1, 4 );
            istringstream atm_line( line.substr( 93, 62 ) );
            atm_line >> atm(0,3)>> atm(0,0) >> atm(0,1) >> atm(0,2) ;                      
            
            // scale from [s] to [m]
            atm(0,1) *= ivg::c;
            atm(0,2) *= ivg::c;           
            atm(0,3) *= ivg::c; 
            // read station and source name
            string ivs_name = remove_spaces_end(line.substr( 48, 8 ));
            string src_name = remove_spaces_end(line.substr( 12, 8 ));
            
            // get epoch
            int y = std::stoi( line.substr( 25,4 ) );
            int m = std::stoi( line.substr( 30,2 ) );
            int d = std::stoi( line.substr( 33,2 ) );
            int h = std::stoi( line.substr( 36,2 ) );
            int min = std::stoi( line.substr( 39,2 ) );
            double s = std::stof( line.substr( 42,4 ) );
            ivg::Date epoch( y,m,d,h,min,s );
            
            // fill raytracing map
            rtmap[ ivs_name ][ epoch ][ src_name ] = atm;    
        }
    }
    inStream.close();
    
    for( auto &sta_it: *trf_ptr )
        sta_it.set_raytracing_data( rtmap[ sta_it.get_name(ivg::staname::ivs_name) ] );
    
#if DEBUG_REFFRAME >=2
   cerr << "--- void raytraced_delays(ivg::Trf * ,  const string )" 
        << " : " << tim.toc() << " s " << endl;
#endif
}

// ...........................................................................
ivg::Matrix c04(const string path, ivg::Date start, ivg::Date end)
{
// ...........................................................................
#ifdef DEBUG_REFFRAME
   cerr << "+++ ivg::Matrix c04(const string , ivg::Date , ivg::Date )" << endl; 
   tictoc tim;
   tim.tic();
#endif
   
    ivg::Matrix data( 1, 11, 0.0 );

//    ifstream inStream;
//    string line,next_line;

    double lod, sigLOD;
    int line_cnt=0;
    
    ifstream inStream(path.c_str(), ios::in);
    if( !inStream.is_open() ){
        throw runtime_error( "ivg::Matrix c04(const string path, ivg::Date start, ivg::Date end): Failed to open file: "+path );
    }else{

        string line;
        while (getline(inStream,line,'\n'))
        {
            // ignore the first 14 header lines of the c04 files
            if(line_cnt > 14)
            {
                double mjd = s2d(line.substr( 13, 6 ));

                if(mjd >= start.get_double_mjd() && mjd <= end.get_double_mjd())
                {
                    //Reading information
                    ivg::Matrix c04( 1, 11 );
                    istringstream c04_line( line.substr( 14, 142 ) );
                    c04_line >> c04(0,0) >> c04(0,1) >> c04(0,2)
                             >> c04(0,3) >> lod >> c04(0,4) >> c04(0,5) >> c04(0,6)
                             >> c04(0,7) >> c04(0,8) >> sigLOD >> c04(0,9) >> c04(0,10) ;
		    
                    c04( 0,3 ) -= ivg::Date(mjd).get_leap_sec();

                    c04( 0,1 ) *= ivg::as2rad; //x pole
                    c04( 0,2 ) *= ivg::as2rad; //y pole
                    c04( 0,3 ) *= ivg::s2rad; //ut1

                    c04( 0,4 ) *= ivg::as2rad; //nut x
                    c04( 0,5 ) *= ivg::as2rad; //nut y

                    c04( 0,6 ) *= ivg::as2rad; //x pole
                    c04( 0,7 ) *= ivg::as2rad; //y pole
                    c04( 0,8 ) *= ivg::s2rad; //ut1

                    c04( 0,9 ) *= ivg::as2rad; //nut x
                    c04( 0,10 ) *= ivg::as2rad; //nut y

                    data.append_rows(c04);
                }
                else if(mjd > end.get_double_mjd())
                    break;
            }
            line_cnt++;
        }
        data.rem_r(0);
    }
    inStream.close();
    
#ifdef DEBUG_REFFRAME
   cerr << "--- ivg::Matrix c04(const string , ivg::Date , ivg::Date )" << " : " << tim.toc() << " s " << endl;
#endif
    return(data);
}

// ...........................................................................
ivg::Matrix cs_erp(const string path, ivg::Date start, ivg::Date end)
{
// ...........................................................................
#ifdef DEBUG_REFFRAME
   cerr << "+++ ivg::Matrix cs_erp(const string , ivg::Date , ivg::Date )" << endl; 
   tictoc tim;
   tim.tic();
#endif
   
    ivg::Matrix data( 1, 11, 0.0 );

    ifstream inStream;
    string line,next_line;

    int i=0;
    while(ivg::parser::get_line(path,inStream, line))
    {
        if(line.substr(0,3) != "EOP" && line.substr(0,1) != "#" )
        {
            double jd = s2d(line.substr( 0, 9 ));
            if(jd >= start.get_jd() && jd <= end.get_jd())
            {
                //Reading information
                ivg::Matrix erp( 1, 11 );
                istringstream erp_line( line.substr( 0, 58 ) );
                // JD, XPO, YPO, UT1, XPO_STD, YPO_STD, UT1_STD
                erp_line >> erp(i,0) >> erp(i,1) >> erp(i,2) >> erp(i,3) >> erp(i,6) >> erp(i,7) >> erp(i,8);

                erp(i,0) -= 2400000.5; // JD to MJD

                erp( i,1 ) *= ivg::as2rad / 10; //x pole
                erp( i,2 ) *= ivg::as2rad / 10; //y pole
                erp( i,3 ) *= ivg::s2rad * 1e-6; //ut1

                erp( i,4 ) = 0.0; //nut x
                erp( i,5 ) = 0.0; //nut y

                erp( i,6 ) *= ivg::as2rad / 10; //x pole std
                erp( i,7 ) *= ivg::as2rad / 10; //y pole std
                erp( i,8 ) *= ivg::s2rad * 1e-6; //ut1 std

                erp( i,9 ) = 0.0; //nut x
                erp( i,10 ) = 0.0; //nut y

                data.append_rows(erp);
            }
            else if(jd > end.get_jd())
                break;
        }
    }
    data.rem_r(0);
    
    inStream.close();
    
    log<DETAIL>("*** ivg::parser::cs_erp Matrix #rows:") % data.rows();
    
#ifdef DEBUG_REFFRAME
   cerr << "--- ivg::Matrix cs_erp(const string , ivg::Date , ivg::Date )" << " : " << tim.toc() << " s " << endl;
#endif
    return(data);
}
// ...........................................................................
ivg::Matrix eops(const string path, ivg::Date start, ivg::Date end)
{
// ...........................................................................
#ifdef DEBUG_REFFRAME
   cerr << "+++ ivg::Matrix eops(const string , ivg::Date , ivg::Date )" << endl; 
   tictoc tim;
   tim.tic();
#endif
   
    ivg::Matrix data( 1, 11, 0.0 );

    ifstream inStream;
    string line,next_line;

    double nutpsi, nuteps;
    int i=0;
    while(ivg::parser::get_line(path,inStream, line))
    {
        if(line.substr(0,3) != "   " && line.substr(0,1) != "#" )
        {
            double mjd = s2d(line.substr( 0, 12 ));
            if(mjd >= start.get_double_mjd() && mjd <= end.get_double_mjd())
            {
                //Reading information
                ivg::Matrix eop( 1, 11 );
                //istringstream eop_line( line.substr( 0, 88 ) );
                istringstream eop_line( line );
                
                // MJD, XPO, YPO, UT1-UTC, NUT_PSI, NUT_EPS, XPO_STD, YPO_STD, UT1-UTC_STD
                eop_line >> eop(i,0) >> eop(i,1) >> eop(i,2) >> eop(i,3) >> nutpsi >> nuteps >> eop(i,6) >> eop(i,7) >> eop(i,8);

                eop( i,1 ) *= ivg::as2rad; //x pole
                eop( i,2 ) *= ivg::as2rad; //y pole
                
                eop( i,3 ) -= ivg::Date(mjd).get_leap_sec(); //ut1
                eop( i,3 ) *= ivg::s2rad; //ut1

                eop( i,4 ) = 0.0; //nut x
                eop( i,5 ) = 0.0; //nut y

                eop( i,6 ) *= ivg::as2rad; //x pole std
                eop( i,7 ) *= ivg::as2rad; //y pole std
                eop( i,8 ) *= ivg::s2rad; //ut1 std

                eop( i,9 ) = 0.0; //nut x std
                eop( i,10 ) = 0.0; //nut y std

/*                if( eop(i,0) == data( data.rows()-1,0 ) )
                {
                   int n = data.rows()-1;
                   for( int m=1; m<4; ++m )
                   {
                      data( n,m ) = ( pow( 1.0/data( n,m+5 ),2.0 )*data( n,m )
                                     +pow( 1.0/eop( 0,m+5 ),2.0 )*eop( 0,m ) )
                                    /( pow( 1.0/data( n,m+5 ),2.0 )+pow( 1.0/eop( 0,m+5 ),2.0 ) );
                      data( n,m+5 ) = 1.0/sqrt(  pow( 1.0/data( n,m+5 ),2.0 )+pow( 1.0/eop( 0,m+5 ),2.0 ) );
                   }
                }
                else*/
                if(eop(0,0)>data(data.rows()-1,0))
                    data.append_rows(eop);
            }
            else if(mjd > end.get_double_mjd())
                break;
        }
    }
    data.rem_r(0);
    
    inStream.close();
    
    log<DETAIL>("*** ivg::parser::eops Matrix #rows:") % data.rows();
    
#ifdef DEBUG_REFFRAME
   cerr << "--- ivg::Matrix eops(const string , ivg::Date , ivg::Date )" << " : " << tim.toc() << " s " << endl;
#endif
    return(data);
}

// ...........................................................................
ivg::Matrix finals(const string path, ivg::Date start, ivg::Date end)
{
// ...........................................................................
#ifdef DEBUG_REFFRAME
   cerr << "+++ ivg::Matrix finals(const string , ivg::Date , ivg::Date )" << endl; 
   tictoc tim;
   tim.tic();
#endif
    
    ivg::Matrix data( 1, 13, 0.0 );

    ifstream inStream;
    string line,next_line;

    double lod, std_lod;
    int i=0;
    while(ivg::parser::get_line(path,inStream, line))
    {
        double mjd = s2d(line.substr( 7, 8 ));

        if(mjd >= start.get_double_mjd() && mjd <= end.get_double_mjd())
        {
            //Reading information
            ivg::Matrix finals( 1, 13 );
            
            if(line.substr(95,1) == " " )
                break;
                    
            finals(i,0) = mjd; // MJD
            finals(i,1) = s2d(line.substr( 18, 9 )) * ivg::as2rad;    // PM-x
            finals(i,6) = s2d(line.substr( 28, 8 )) * ivg::as2rad;   // PM-x error
            
            finals(i,2) = s2d(line.substr( 37, 9 )) * ivg::as2rad;   // PM-y
            finals(i,7) = s2d(line.substr( 47,8 )) * ivg::as2rad;   // PM-y error
            
            finals(i,3) = ( s2d(line.substr( 58, 10 )) - ivg::Date(mjd).get_leap_sec() ) * ivg::s2rad;   // UT1-UTC
            finals(i,8) = s2d(line.substr( 69,9 )) * ivg::s2rad;  // UT1-UTC error
            
//            finals(i,11) = ( s2d(line.substr( 79,9 )) / 1000 ) * ivg::s2rad;   // LOD
//            finals(i,12) = ( s2d(line.substr( 87,7 )) / 1000 ) * ivg::s2rad;   // LOD error
            
            finals(i,4) = ( s2d(line.substr( 98,8 )) / 1000 ) * ivg::as2rad;   // dX
            finals(i,9) = ( s2d(line.substr( 108,7 )) / 1000 ) * ivg::as2rad; // dX error

            finals(i,5) = ( s2d(line.substr( 117,8 )) / 1000 ) * ivg::as2rad; // dY
            finals(i,10) = ( s2d(line.substr( 127,7 )) / 1000 ) * ivg::as2rad; // dY error

            data.append_rows(finals);

        }
        else if(mjd > end.get_double_mjd())
            break;
    }
    data.rem_r(0);
    
    inStream.close();
    
#ifdef DEBUG_REFFRAME
   cerr << "--- ivg::Matrix finals(const string , ivg::Date , ivg::Date )" << " : " << tim.toc() << " s " << endl;
#endif
    return(data);
}

// ...........................................................................
ivg::Matrix igs_erp(const string path, ivg::Date start, ivg::Date end)
// ...........................................................................
{
#ifdef DEBUG_REFFRAME
   cerr << "+++ ivg::Matrix igs_erp(const string , ivg::Date , ivg::Date )" << endl; 
   tictoc tim;
   tim.tic();
#endif    
    ivg::Matrix data( 1, 13, 0.0 );

    ifstream inStream;
    string line;

    int i=0;
    while(ivg::parser::get_line(path,inStream, line))
    {
        if(line.substr(0,14) == "EOP  SOLUTION" )
            break;
    }
    while(ivg::parser::get_line(path,inStream, line))
    {
        if(line.substr(0,2) != "  " )
        {
            double mjd = s2d(line.substr( 0, 8 ));

            if(mjd >= start.get_double_mjd() && mjd <= end.get_double_mjd())
            {
                //Reading information
                ivg::Matrix erp( 1,13,0.0 );

                erp(i,0) = mjd; // MJD
                erp(i,1) = s2d(line.substr( 10, 8 ))*1e-6*ivg::as2rad;    // PM-x
                erp(i,6) = s2d(line.substr( 44, 6 ))*1e-6*ivg::as2rad;   // PM-x error

                erp(i,2) = s2d(line.substr( 18, 8 ))*1e-6*ivg::as2rad;   // PM-y
                erp(i,7) = s2d(line.substr( 51,6 ))*1e-6*ivg::as2rad;   // PM-y error

                erp(i,3) = ( s2d(line.substr( 27, 8 ))*1e-6-ivg::Date(mjd).get_leap_sec() )*ivg::s2rad;   // UT1-UTC
                erp(i,8) = s2d(line.substr( 58,6 ))*1e-6*ivg::s2rad;  // UT1-UTC error

    //            erp(i,11) = 0.0 * ivg::s2rad; // LOD
    //            erp(i,12) = 0.0 * ivg::s2rad; // LOD error
    //            
    //            erp(i,4) = 0.0 * ivg::s2rad;  // dX
    //            erp(i,9) = 0.0 * ivg::s2rad;  // dX error
    //
    //            erp(i, 5) = 0.0 * ivg::s2rad; // dY
    //            erp(i,10) = 0.0 * ivg::s2rad; // dY error

                data.append_rows(erp);

            }
            else if(mjd > end.get_double_mjd())
                break;
        }
    }
    data.rem_r(0);
    
    inStream.close();
    
#ifdef DEBUG_REFFRAME
   cerr << "--- ivg::Matrix igs_erp(const string , ivg::Date , ivg::Date )" << " : " << tim.toc() << " s " << endl;
#endif
    return(data);
}
// ...........................................................................
map<string,string> correspondence(const string path)
// ...........................................................................
{
#if DEBUG_REFFRAME >=2
   cerr << "+++ map<string,string> correspondence(const string)" << endl; 
   tictoc tim;
   tim.tic();
#endif
   
    map<string,string> corres;

    ifstream inStream;
    string line;

    while( ivg::parser::get_line(path, inStream, line))
    {
        if(line.substr(0,4) == "VLBI")
            corres[remove_spaces_end(line.substr( 6, 8 ))] = remove_spaces_end(line.substr( 18, 8 ));
    }
    
    inStream.close();
    
#if DEBUG_REFFRAME >=2
    cerr << "--- map<string,string> correspondence(const string)" << " : " << tim.toc() << " s " << endl;
#endif 
    return(corres);
}

// ...........................................................................
vector< map<ivg::staname,string> >  nscodes_parser(
    const vector<string> station_names, const ivg::staname type,
    const string path)
// ...........................................................................
{
#if DEBUG_REFFRAME >=2
   cerr << "+++ vector< map<ivg::staname,string> >  nscodes_parser(const vector<string> , const ivg::staname ,const string , map<string,string> )" << endl; 
   tictoc tim;
   tim.tic();
#endif
   
    vector<map<ivg::staname,string>> ns_codes;

    ifstream inStream;
    string line;

    ivg::parser::get_line(path, inStream, line);

    while( ivg::parser::get_line(path, inStream, line) )
    {
        if(line.substr(0,1) != "*")
        {
            map<ivg::staname,string> names;
            names[ivg::staname::lettercode] = line.substr( 1,2 );
            names[ivg::staname::ivs_name] = remove_spaces_end(line.substr( 4,8 ));
            names[ivg::staname::domes_no] = line.substr( 13,9 );
            names[ivg::staname::cdp] = line.substr( 23,4 );
            names[ivg::staname::description] = line.substr( 28 );
//            names[ivg::staname::corres] = correspondence[names[ivg::staname::ivs_name]];

            // get ALL stations in case of empty vector
            if(station_names.empty())
                ns_codes.push_back(names);
            else
            {
                for(auto &sta: station_names)
                    if(sta == names[type])
                        ns_codes.push_back(names);
            }
        }
    }
    
    inStream.close();
    
#if DEBUG_REFFRAME >=2
    cerr << "--- vector< map<ivg::staname,string> >  nscodes_parser(const vector<string> , const ivg::staname ,const string , map<string,string> )" << " : " << tim.toc() << " s " << endl;
#endif 
    return(ns_codes);
}

// ...........................................................................
vector<ivg::Analysis_station> ssc_parser(const string path, const vector< map<ivg::staname,string> > &nscodes)
// ...........................................................................
{
#if DEBUG_REFFRAME >=2
   cerr << "+++ vector<ivg::Analysis_station> ssc_parser(const string ,  const vector< map<ivg::staname,string> > )" << endl; 
   tictoc tim;
   tim.tic();
#endif
   
    // _stations vector containing all station-specific information
    vector<ivg::Analysis_station> stations;

    ifstream inStream;
    string line;
    ivg::parser::get_line(path, inStream, line);

    // first line contains overall reference epoch in format 2008.0
    int eppos=line.find("STATION POSITIONS AT EPOCH")+27;
    ivg::Date _ref_epoch( atoi( line.substr( eppos,6 ).c_str() ), 1.0 );
    if(_ref_epoch.get_double_mjd() < 0)
        throw runtime_error( "vector<ivg::Analysis_station> ssc_parser(...): Unexpected reference epoch in SSC-File: "+path );

    while( ivg::parser::get_line(path, inStream, line) )
    {
        int org_linesize=line.size();
        line=remove_spaces_end(line);
        // lines containing data start with int
        istringstream inpStream( line.substr( 0, 1 ) );
        int inpValue;
        if (inpStream >> inpValue)
        {
            string domes_no = line.substr( 0,9 );
	    
            // search domes_no in existing nscodes
            for (int i=0 ; i<nscodes.size(); i++)
            {
	        
                if(nscodes.at(i).at(ivg::staname::domes_no) == domes_no)
                {       
                    bool new_station;
                    string longline = "";
                    
                    if( line.size() > 100  )
                        longline = line;
                    
                    int access = 0;
                    
                    if(line.substr( 85, 1 ) == "." && line.size() > 100 && line.substr(101,1) != "1")
                        access = 103;
                    else if(line.substr( 85, 1 ) == "." && line.size() > 100 && line.substr(101,1) == "1")
                        access = -1;
                    else if(line.substr( 85, 1 ) == "." && line.size() <= 100 )
                        access = -1;
                    else if(line.substr( 85, 1 ) != "." && line.size() > 94)
                        access = 97;
                    else if(line.substr( 85, 1 ) != "." && line.size() <= 94)
                        access = -1;

		    
                    ivg::Matrix pos0( 3, 1 );
                    istringstream pos_line( line.substr( 37, 38 ) );
                    pos_line >> pos0( 0 ) >> pos0( 1 ) >> pos0( 2 );

                    ivg::parser::get_line(path, inStream, line);
                    Matrix vel0( 3, 1 );
                    istringstream vel_line( line.substr( 37, 38 ) );
                    vel_line >> vel0( 0 ) >> vel0( 1 ) >> vel0( 2 );
                    
                    if( access == -1 )
                    {
                        vector<ivg::Date> discont = {ivg::Date ( "70:001:00000" , "SINEX" )};
                        ivg::Analysis_station new_station(pos0, vel0, _ref_epoch, discont, nscodes.at(i));
                        stations.push_back(new_station);
                    }
                    else
                    {
                        vector<ivg::Analysis_station>::iterator iter;
                        for(iter=stations.begin(); iter < stations.end(); iter++)
                        {
                            if((*iter).get_name(ivg::staname::domes_no) == domes_no)
                               (*iter).add_discontinuity(pos0, vel0, _ref_epoch, ivg::Date ( longline.substr(access,12) , "SINEX" ));
                        }
                    }
                    
                    break;
                } 
            }
        }
    }
    
    inStream.close();
        
#if DEBUG_REFFRAME >=2
    cerr << "--- vector<ivg::Analysis_station> ssc_parser(const string ,  const vector< map<ivg::staname,string> > )" << " : " << tim.toc() << " s " << endl;
#endif 
    return(stations);
}
// ...........................................................................
void stations_cat(ivg::Trf * trf_ptr,  const string path)
{
// ...........................................................................
#if DEBUG_REFFRAME >=2
   cerr << "+++ void stations_cat(ivg::Trf * trf_ptr,  const string path)" << endl; 
   tictoc tim;
   tim.tic();
#endif
   
    // only setting the ivg::staname::ant_name in this function
    ifstream inStream;
    string line;

    while( ivg::parser::get_line(path, inStream, line) )
    {
        if(line.substr(0,1) != "*")
        {
            string sked_antenna_name = remove_spaces_end(line.substr( 6, 8 ));
            string ivs_name = remove_spaces_end(line.substr( 21, 8 ));
//            string rack_type = remove_spaces_end(line.substr( 38, 4 )); // NOT USED RIGHT NOW
            
            ivg::Analysis_station * station;
            if(trf_ptr->get_station(&station, ivs_name, ivg::staname::lettercode))
                station->set_name(ivg::staname::ant_name, sked_antenna_name);
      }
   }
    
   inStream.close();
   
#if DEBUG_REFFRAME >=2
    cerr << "--- void stations_cat(ivg::Trf * trf_ptr,  const string path)" << " : " << tim.toc() << " s " << endl;
#endif 
}
// ...........................................................................
void antenna_cat(ivg::Trf * trf_ptr,  const string path)
{
// ...........................................................................
#if DEBUG_REFFRAME >=2
   cerr << "+++ void antenna_cat(ivg::Trf * trf_ptr,  const string path)" << endl; 
   tictoc tim;
   tim.tic();
#endif
   
    ifstream inStream;
    string line;

    char char_begin = 'A';
    vector<string> used_ids;
    while( ivg::parser::get_line(path, inStream, line) )
    {
        if(line.substr(0,1) != "*")
        {
            // tokenize each row
            vector<string> tokens = get_tokens(line);
            if( tokens.size() > 1 )
            {
                string ant_name = tokens.at( 1 );
                
                // in case of TIGO we need some special handling
                if( ant_name == "TIGO" )
                    ant_name = "TIGOCONC";
                
                string equip_id = tokens.at( 14 );
                
                ivg::Analysis_station * station;
                // changed from ivg::staname::equip_id to MINSTA
                // should be no problem, because corres not set anymore
                if(trf_ptr->get_station(&station, ant_name, ivg::staname::MINSTA, ivg::staname::MAXSTA))
                {
                    // some ids are used multiple times (see antenna.cat)
                    // therefore we need to determine a new character
                    string new_id = tokens.at( 0 );
                    if(find( used_ids.begin(), used_ids.end(), new_id) == used_ids.end()) 
                    {
                        station->set_antenna_info().id = new_id;
                        used_ids.push_back(new_id);
                    }
                    else
                    {
                        while(1)
                        {
                            stringstream ss;
                            ss << char_begin;
                            ss >> new_id;
                            if(find( used_ids.begin(), used_ids.end(), new_id) == used_ids.end())
                            {
                                used_ids.push_back(new_id);
                                break;
                            }
                            else                            
                                char_begin++;  // -> A -> B -> C -> D 
                        }
                        station->set_antenna_info().id = new_id;
                    }
                    station->set_antenna_info().axis_type = tokens.at( 2 );
                    station->set_antenna_info().num_id = tokens.at( 14 );
                    station->set_antenna_info().offset = s2d( tokens.at( 3 ) );
                    station->set_antenna_info().ant_diameter = s2d( tokens.at( 12 ) );
                    station->set_antenna_info().azi_rate =  s2d( tokens.at( 4 ) )*(M_PI/180.0)/60.0; // from degree/min to rad/sec 
                    station->set_antenna_info().ele_rate =  s2d( tokens.at( 8 ) )*(M_PI/180.0)/60.0; // from degree/min to rad/sec
                    station->set_antenna_info().azi_const =  s2d( tokens.at( 5 ) ); // sec
                    station->set_antenna_info().ele_const =  s2d( tokens.at( 9 ) ); // sec
                    station->set_antenna_info().azi_min =  s2d( tokens.at( 6 ) )*(M_PI/180.0); // rad
                    station->set_antenna_info().azi_max =  s2d( tokens.at( 7 ) )*(M_PI/180.0); // rad
                    station->set_antenna_info().ele_min =  s2d( tokens.at( 10 ) )*(M_PI/180.0); // rad
                    station->set_antenna_info().ele_max =  s2d( tokens.at( 11 ) )*(M_PI/180.0);  // rad
                    
                    // if first limit is negative, correct it for 360° (correspondes to sked procedure)
                    if(station->set_antenna_info().azi_min < 0.0)
                    {
                        station->set_antenna_info().azi_min += 2*M_PI;
                        station->set_antenna_info().azi_max += 2*M_PI;
                    }
                    // setting A_line for easy sked-file writing
                    station->set_antenna_info().A_line = line.substr(17);
                    station->set_name(ivg::staname::equip_id, equip_id);
                }
            }
        }
    }
    
    inStream.close();

#if DEBUG_REFFRAME >=2
    cerr << "--- void antenna_cat(ivg::Trf * trf_ptr,  const string path)" << " : " << tim.toc() << " s " << endl;
#endif 
}
// ...........................................................................
void flux_cat(ivg::Crf * crf_ptr,  const string path, bool isSkedFile)
{
// ...........................................................................
#if DEBUG_REFFRAME >=2
   cerr << "+++ void flux_cat(ivg::Crf * crf_ptr,  const string path,  bool isSkedFile)" << endl; 
   tictoc tim;
   tim.tic();
#endif
   
    ifstream inStream;
    string line;

    if(isSkedFile){
        while( ivg::parser::get_line(path, inStream, line) )
        {
            if ( boost::algorithm::starts_with(line, "$FLUX")) {
                while (ivg::parser::get_line(path, inStream, line) && !boost::algorithm::starts_with(line, "$")) {
                    ivg::parser::flux_cat_line(crf_ptr, line);
                }
                break;
            }
        }      
    } else {
        while( ivg::parser::get_line(path, inStream, line) )
        {
            flux_cat_line(crf_ptr, line);
        }
    }

    
    inStream.close();
    
#if DEBUG_REFFRAME >=2
    cerr << "--- void flux_cat(ivg::Crf * crf_ptr,  const string path)" << " : " << tim.toc() << " s " << endl;
#endif 
}

// ...........................................................................
void flux_cat_line(ivg::Crf * crf_ptr,  std::string line)
{
// ...........................................................................
#if DEBUG_REFFRAME >=2
   cerr << "+++ void flux_cat_line(ivg::Crf * crf_ptr,  std::string line)" << endl; 
   tictoc tim;
   tim.tic();
#endif
    if(line.substr(0,1) != "*")
    {
        // tokenize each row
        vector<string> tokens = get_tokens(line);
        if( tokens.size() > 1 )
        {
           string name = tokens.at( 0 ); // NRAO140
           string band   = tokens.at( 1 ); // X or S
           string type = tokens.at( 2 ); // M or B

           ivg::Source * source;
           if(crf_ptr->get_source(&source, name))
           {
                ivg::band current_band;
                ivg::flux_info current_flux_info;
                current_flux_info.flux_line = line; // save complete line for $FLUX block in skd file
                if(band == "X")
                    current_band = ivg::band::X;
                else if(band == "S")
                    current_band = ivg::band::S;

                current_flux_info.type = type;  // M or B
                //* Source Band Type Flux   MajAx  Ratio  PA   Off1   Off2
                //* Name        M    (Jy)   (mas)              (mas)  (mas)
                //0003-066  X   M    3.82   0.80   0.40   10.   0.0   0.0
                if(type == "M" && tokens.size() == 9)
                {
                    current_flux_info.flux = s2d(tokens.at( 3 ));
                    current_flux_info.major_axis = s2d(tokens.at( 4 ));
                    current_flux_info.ratio = s2d(tokens.at( 5 ));
                    current_flux_info.pa = s2d(tokens.at( 6 ));
                    current_flux_info.off1 = s2d(tokens.at( 7 ));
                    current_flux_info.off2 = s2d(tokens.at( 8 ));
                }
                //* Source Band Type 0.0  Flux   Baseline Flux   Baseline
                //* Name          B       (Jy)    limit    (Jy)   limit
                //   3C84   X     B  0.0  10.00   600.0    0.00   13000.0
                else if(type == "B")
                {
                    for( int i=4; i<tokens.size(); i+=2)
                    {
                        ivg::Matrix new_row = ivg::Matrix(vector<double>({s2d(tokens.at(i+1)), s2d(tokens.at(i))})).transpose();
                        current_flux_info.baseline_flux.append_rows(new_row);
                    }
                }
               source->add_band_flux_info(current_band, current_flux_info);
           }
        }
    }
#if DEBUG_REFFRAME >=2
    cerr << "--- void flux_cat_line(ivg::Crf * crf_ptr,  std::string line)" << " : " << tim.toc() << " s " << endl;
#endif 
}


// ...........................................................................
void equip_cat(ivg::Trf * trf_ptr,  const string path)
{
// ...........................................................................
#if DEBUG_REFFRAME >=2
   cerr << "+++ void equip_cat(ivg::Trf * trf_ptr,  const string path)" << endl; 
   tictoc tim;
   tim.tic();
#endif
   
    ifstream inStream;
    string line;
    while( ivg::parser::get_line(path, inStream, line) )
    {
        if(line.substr(0,1) != "*")
        {
            // tokenize each row            
            vector<string> tokens = get_tokens(line);
            if( tokens.size() > 1 )
            {
                string ant_name = tokens.at( 0 );
                string equip_id = tokens.at( 1 );
                
                ivg::Analysis_station * station;                           
                if(trf_ptr->get_station(&station, equip_id, ivg::staname::equip_id, ivg::staname::MAXSTA))
                {
                    if( tokens.size() > 6 && tokens.at(5) != "C" && tokens.at(7) != "C")
                    {
                       ivg::Equip equip_info;
                       equip_info.head = tokens.at( 3 );
                       equip_info.tape_length = tokens.at( 4 );
                       vector<double> x_sefd = {s2d( tokens.at( 6 ) )};
                       vector<double> s_sefd = {s2d( tokens.at( 8 ) )};

                       if( tokens.size() > 9 && tokens.at( 9 ) != "" )
                       {
                           equip_info.recorder = tokens.at(tokens.size()-1);
                           
                            // get SEFD parameter [X1,X2,X3,S1,S2,S3]
                            if( tokens.at( 9 ) == "X" )
                            {
                                x_sefd.push_back(s2d( tokens.at( 10 )));
                                x_sefd.push_back(s2d( tokens.at( 11 )));
                                x_sefd.push_back(s2d( tokens.at( 12 )));
                                s_sefd.push_back(s2d( tokens.at( 14 )));
                                s_sefd.push_back(s2d( tokens.at( 15 )));
                                s_sefd.push_back(s2d( tokens.at( 16 )));
                            }
                            else if( tokens.at( 9 ) == "S" )
                            {
                                s_sefd.push_back(s2d( tokens.at( 10 )));
                                s_sefd.push_back(s2d( tokens.at( 11 )));
                                s_sefd.push_back(s2d( tokens.at( 12 )));
                                x_sefd.push_back(s2d( tokens.at( 14 )));
                                x_sefd.push_back(s2d( tokens.at( 15 )));
                                x_sefd.push_back(s2d( tokens.at( 16 )));
                            }                                        
                       }
                       
                       // e.g. " 33  WETTZELL  2x56000 17640   X   750   S  1115 S 1.0 0.934 0.0660 X 1.0 0.948 0.0516 Mark4 MARK5A"
                       equip_info.line = line.substr(11); // for easy skd-file writing
                       
                       map< ivg::band, vector<double> > sefd;
                       sefd[ivg::band::X] = x_sefd;
                       sefd[ivg::band::S] = s_sefd;
                       
                       equip_info.sefd = sefd;
                       station->add_equip_info(equip_info);
                    }
                }
            }
        }
    }
    
    inStream.close();

#if DEBUG_REFFRAME >=2
    cerr << "--- void equip_cat(ivg::Trf * trf_ptr,  const string path)" << " : " << tim.toc() << " s " << endl;
#endif 
}
// ...........................................................................
string freq_cat(ivg::Trf * trf_ptr,  const string path, string freqname)
{
// ...........................................................................
#if DEBUG_REFFRAME >=2
   cerr << "+++ string freq_cat(ivg::Trf * trf_ptr,  const string path)" << endl; 
   tictoc tim;
   tim.tic();
#endif
   
    string rxname; // e.g. SX_WIDE
    
    ifstream inStream;
    string line;
    while( ivg::parser::get_line(path, inStream, line) )
    {
        if(line.substr(0,1) != "*" && line.substr(0,1) != "!" && line.size() > 10)
        {
            vector<string> tokens = get_tokens(line);
            if(tokens.at(0)== freqname)
            {
                rxname = tokens.at(3);
                string freq_lc = tokens.at(1); // frequency lettercode, e.g. 8F
                ivg::Matrix x_freqs, s_freqs;
                map<int, ivg::Frequency> tmp_info;
                while( ivg::parser::get_line(path, inStream, line) && line.substr(0,1) == "-" )
                {
                    double value = s2d(line.substr(11,7));
                    if(line.substr(2,1) == "X")
                        x_freqs.append_rows(value);
                    else if(line.substr(2,1) == "S")
                        s_freqs.append_rows(value);
                    
                    tokens = get_tokens(line);
                    int channel_id = stoi(tokens.at(5).substr(2));
                    int bbc = stoi(tokens.at(6));
                    for( vector<ivg::Analysis_station>::iterator sta_it = trf_ptr->begin(); sta_it != trf_ptr->end(); sta_it++ ) 
                    {
                        sta_it->set_channel_setup().freq_sq_name = freqname;
                        sta_it->set_channel_setup().freq_lc = freq_lc;
                        sta_it->set_channel_setup().cha_bbc[channel_id] = bbc;
                        sta_it->set_channel_setup().freq_info[bbc].band = tokens.at(1);
                        sta_it->set_channel_setup().freq_info[bbc].pol = tokens.at(2);
                        sta_it->set_channel_setup().freq_info[bbc].frequency = s2d(tokens.at(3));
                        sta_it->set_channel_setup().freq_info[bbc].sideband = tokens.at(4);
                        sta_it->set_channel_setup().freq_info[bbc].channel = tokens.at(5);
                        sta_it->set_channel_setup().freq_info[bbc].phase_cal_freq = s2d(tokens.at(7));
                    }

                } 

                // right now we save this information redundantly
                for(auto &station: (*trf_ptr))
                {
                    station.add_band_frequency_sequence(ivg::band::X, x_freqs);
                    station.add_band_frequency_sequence(ivg::band::S, s_freqs);
                }

                break;
            }
        }
    }
    
    inStream.close();

    return rxname;
#if DEBUG_REFFRAME >=2
    cerr << "--- string freq_cat(ivg::Trf * trf_ptr,  const string path)" << " : " << tim.toc() << " s " << endl;
#endif 
}
// ...........................................................................
void loif_cat(ivg::Trf * trf_ptr,  const string path, const string rx_path, string rx_name)
{
// ...........................................................................
#if DEBUG_REFFRAME >=2
   cerr << "+++ void loif_cat(ivg::Trf * trf_ptr,  const string path, const string rx_path, string rx_name)" << endl; 
   tictoc tim;
   tim.tic();
#endif
   
    // in order to be able to load loif_cat we need the relation from rx-file
    // e.g. SX_WIDE -> CDP_WIDE for WETTZELL
    ifstream inStream;
    string line;
    map<ivg::Analysis_station*, string> loif_name;
    while( ivg::parser::get_line(rx_path, inStream, line) )
    {
        if(line.substr(0,1) != "*")
        {         
            vector<string> tokens = get_tokens(line);
            if( tokens.size() >= 1 && rx_name == tokens.at(0))
            {
                while( ivg::parser::get_line(rx_path, inStream, line) && ( line.substr(1,1) == "-" || line.substr(2,1) == "-" ) )
                {
                    tokens = get_tokens(line);
                    ivg::Analysis_station * station;                           
                    if(tokens.size() >= 2 && trf_ptr->get_station(&station, tokens.at(1), ivg::staname::MINSTA, ivg::staname::MAXSTA))
                        loif_name[station] = tokens.at(2);
                }
                break;
            }
        }
    }
    inStream.close(); // from rx.cat file
    
    if(loif_name.size() != trf_ptr->get_number_stations())
        throw runtime_error("ERROR: "+rx_path+" does not contain information about all stations to be scheduled.");
    
    ifstream inStream_2;
    // opens file as many times as stations; but we have a break!
    for(auto &loif: loif_name)
    {
        // now loif_name contains e.g. (pseudocode) loif_name[WETTZELL] = CDP_WIDE in case of rx_name == SX_WIDE
        while( ivg::parser::get_line(path, inStream_2, line) )
        {
            if(line.substr(0,1) != "*")
            {
                vector<string> tokens = get_tokens(line);
                if( tokens.at(0) == loif.second)
                {
                    map<int, ivg::Frequency> tmp_info;
                    while( ivg::parser::get_line(path, inStream_2, line) && line.substr(1,1) == "-" )
                    {
                        tokens = get_tokens(line);
                        int bbc = stoi(tokens.at(1));
                        loif.first->set_channel_setup().loif_name = loif.second;
                        loif.first->set_channel_setup().freq_info[bbc].if_channel = tokens.at(2);
                        loif.first->set_channel_setup().freq_info[bbc].lo_frequency = s2d(tokens.at(4));
                        loif.first->set_channel_setup().freq_info[bbc].lo_sideband = tokens.at(5);
                    }
                    
                    break;
                }
            }
        }
        inStream_2.close();  // from loif.cat file
    }
   
#if DEBUG_REFFRAME >=2
    cerr << "--- void loif_cat(ivg::Trf * trf_ptr,  const string path, const string rx_path, string rx_name)" << " : " << tim.toc() << " s " << endl;
#endif 
}
// ...........................................................................
void tracks_cat(ivg::Trf * trf_ptr,  const string path, const string rec_path, string rec_mode)
// ...........................................................................
{
#if DEBUG_REFFRAME >=2
   cerr << "+++ void tracks_cat(ivg::Trf * trf_ptr,  const string path, const string rec_path, string rec_mode)" << endl; 
   tictoc tim;
   tim.tic();
#endif
    
    // in order to be able to load tracks_cat we need the relation from rec-file
    ifstream inStream;
    string line;
    map<ivg::Analysis_station*, string> tracks_name,recs_name,recs_hdpos_names;
    while( ivg::parser::get_line(rec_path, inStream, line) )
    {
        if(line.substr(0,1) != "*")
        {         
            vector<string> tokens = get_tokens(line);
            if( tokens.size() >= 1 && rec_mode == tokens.at(0))
            {
                while( ivg::parser::get_line(rec_path, inStream, line) && line.size()>0 && line.substr(1,1) == "-" )
                {
                    tokens = get_tokens(line);
                    ivg::Analysis_station * station;                           
                    if(tokens.size() >= 2 && trf_ptr->get_station(&station, tokens.at(1), ivg::staname::MINSTA, ivg::staname::MAXSTA))
                    {
                        recs_hdpos_names[station] = tokens.at(2);
                        tracks_name[station] = tokens.at(3);
                        recs_name[station] = tokens.at(4);
                    }
                }
                break;
            }
        }
    }
    inStream.close(); // from rec_path file
    
    if(recs_name.size() != trf_ptr->get_number_stations())
        throw runtime_error("ERROR: "+rec_path+" does not contain information about all stations to be scheduled.");
       
    // now tracks_name contains e.g. (pseudocode) tracks_name[WETTZELL] = 14U2L-2-1 in case of rec_mode = 00-16-0-1
    // furthermore recs_name[WETTZELL] = Mk34
    ifstream inStream_2;
    // opens file as many times as stations; but we have a break!
    for(auto &track: tracks_name)
    {
        while( ivg::parser::get_line(path, inStream_2, line) )
        {
            if(line.substr(0,1) != "*")
            {
                vector<string> tokens = get_tokens(line);
                if( tokens.size() > 0 && tokens.at(0) == track.second)
                {
                    string fanout_fac = tokens.at(2)+":"+tokens.at(1);
                    map<int, ivg::Frequency> tmp_info;
                    while( ivg::parser::get_line(path, inStream_2, line) && line.size()>0 && line.substr(1,1) == "-" )
                    {
                        tokens = get_tokens(line);
                        int channel_id = stoi(tokens.at(1));
                        int bbc = track.first->set_channel_setup().cha_bbc[channel_id];
                        track.first->set_channel_setup().freq_info[bbc].track_info = tokens.at(2);
                        track.first->set_channel_setup().freq_info[bbc].fanout_fac = fanout_fac;
                        track.first->set_channel_setup().rec_format = recs_name[track.first];
                        track.first->set_channel_setup().rec_hdpos_format = recs_hdpos_names[track.first];
                    }
                    break;
                }
            }
        }
        inStream_2.close();  // from tracks.cat file
    }
   
#if DEBUG_REFFRAME >=2
    cerr << "--- void tracks_cat(ivg::Trf * trf_ptr,  const string path, const string rec_path, string rec_mode)" << " : " << tim.toc() << " s " << endl;
#endif 
}
// ...........................................................................
void hdpos_cat(ivg::Trf * trf_ptr,  const string path)
// ...........................................................................
{
#if DEBUG_REFFRAME >=2
   cerr << "+++ void hdpos_cat(ivg::Trf * trf_ptr,  const string path)" << endl; 
   tictoc tim;
   tim.tic();
#endif
    
    // ugly but who care?!
    // stores information for each analysis center, e.g:
    //  MK3V-A
    // -   11(-319) 21(31)   31(-271) 41(79)   51(-223) 61(127)
    // -   71(-175) 81(175)  91(-127) A1(223)  B1(-79)  C1(271)
    // -   D1(-31)  E1(319) 
    for( vector<ivg::Analysis_station>::iterator sta_it = trf_ptr->begin(); sta_it != trf_ptr->end(); sta_it++ ) 
    {
        ifstream inStream;
        string line, hdpos_lines,hdpos_format;
        while( ivg::parser::get_line(path, inStream, line) )
        {
            if(line.substr(0,1) != "*")
            {
                vector<string> tokens = get_tokens(line);
                if(tokens.size()>0)
                {
                    if(sta_it->set_channel_setup().rec_hdpos_format == tokens.at(0))
                    {
                        string hdpos_lines;
                        while( ivg::parser::get_line(path, inStream, line) && line.substr(1,1) == "-" )
                            hdpos_lines += sta_it->get_antenna_info().id + " " + sta_it->set_channel_setup().freq_lc + line.substr(2) + "\n";
                        
                        sta_it->set_channel_setup().hdpos_lines = hdpos_lines.substr(0,hdpos_lines.size()-2); // remove new line
                        break;
                    }
                }
            }
        }
        inStream.close(); 
    }
  
#if DEBUG_REFFRAME >=2
    cerr << "--- void hdpos_cat(ivg::Trf * trf_ptr,  const string path)" << " : " << tim.toc() << " s " << endl;
#endif 
}
// ...........................................................................
void mask_cat(ivg::Trf * trf_ptr,  const string path)
{
// ...........................................................................
#if DEBUG_REFFRAME >=2
   cerr << "+++ void mask_cat(ivg::Trf * trf_ptr,  const string path)" << endl; 
   tictoc tim;
   tim.tic();
#endif
   
    ifstream inStream;
    string line;
    while( ivg::parser::get_line(path, inStream, line) )
    {
        if(line.substr(0,1) != "*" && line.size() > 2)
        {
            string full_line = line;
            vector<string> tokens = get_tokens(line);
            // find station from current line in trf
            ivg::Analysis_station * station;
            if(trf_ptr->get_station(&station, tokens.at(1), ivg::staname::ivs_name, ivg::staname::lettercode))
            {
                // in case of H-mask
                if(tokens.at(0) == "H")
                {
                    // read as many as lines are given
                    while( ivg::parser::get_line(path, inStream, line) && line.substr(0,1) != "*" && line.substr(1,1) != "C" && line.substr(1,1) != "H")
                        full_line += " " + line.substr(2);
                    
                    // split the full_line into azimut and elevation and save it in a map
                    // min_elevation correction and interpolation need to be performed in Analysis_station
                    ivg::Matrix mask;
                    vector<string> azi_ele = get_tokens(full_line);
                    // check if we always have azi-ele pairs, otherwise correct it by closing the 360°! (e.g. BADARY)
                    if((azi_ele.size()-3)%2 != 0)
                        azi_ele.push_back(azi_ele.at(4));
                    
                    for(int i=3; i<azi_ele.size(); i+=2)
                    {
                        ivg::Matrix row(1,2,0.0);
                        row(0,0) = s2d(azi_ele.at(i))*(M_PI/180.0);
                        row(0,1) = s2d(azi_ele.at(i+1))*(M_PI/180.0);
                        mask.append_rows(row);
                    }
                    station->add_mask_info(mask);
                    station->set_antenna_info().H_line = full_line.substr(15); // for easy skd-file writing
                }
                else if(tokens.at(0) == "C")
                {
                    throw runtime_error("void mask_cat(...): mask.cat Type C not supported yet.");
                    // NOT IMPLEMENTED YET
                }
            }
        }
    }
    
    inStream.close();

#if DEBUG_REFFRAME >=2
    cerr << "--- void mask_cat(ivg::Trf * trf_ptr,  const string path)" << " : " << tim.toc() << " s " << endl;
#endif 
}
// ...........................................................................
void sp3(ivg::Crf * crf_ptr,  const string path)
{
// ...........................................................................
#if DEBUG_REFFRAME >=2
   cerr << "+++ void sp3(ivg::Crf * crf_ptr,  const string path)" << endl; 
   tictoc tim;
   tim.tic();
#endif
   
    ifstream inStream;
    string line;
    
    ivg::Date epoch;
    map<ivg::Source *, ivg::Matrix> sp3_data;
    while( ivg::parser::get_line(path,inStream, line))
    {
        if(line.substr(0,1) == "*")
        {
            vector<string> time = get_tokens(line); 
            
            int y = std::stoi( time.at(1) );
            int m = std::stoi( time.at(2) );
            int d = std::stoi( time.at(3) );
            int h = std::stoi( time.at(4) );
            int min = std::stoi( time.at(5) );
            double sec = s2d( time.at(6) );
            
            epoch =  ivg::Date(y,m,d,h,min,sec);
            ivg::parser::get_line(path,inStream, line);
        }
        
        string sat_name = line.substr(0,4); // e.g. PG02
        
        ivg::Matrix new_column(4,1);
        new_column(0) = epoch.get_double_mjd(); // keep epoch in GPS-Time
        new_column(1) = s2d(line.substr(4,14));
        new_column(2) = s2d(line.substr(18,14));
        new_column(3) = s2d(line.substr(32,14));
        
        ivg::Source * source;
        if(crf_ptr->get_source( &source, sat_name ))
            sp3_data[source].append_cols(new_column);
    }
    
    // set final sp3 data containing orbit information for each satellite
    for(auto &data: sp3_data)
        data.first->set_sp3(data.second.transpose());
    
    inStream.close();
   
    log<INFO>("*** ivg::parser::sp3 Loaded sp3-IGS-final-orbits from ") % path;
    
#if DEBUG_REFFRAME >=2
    cerr << "--- void void sp3(ivg::Crf * crf_ptr,  const string path)" << " : " << tim.toc() << " s " << endl;
#endif 
}
// ...........................................................................
void tle(ivg::Crf * crf_ptr,  const string path)
{
// ...........................................................................
#if DEBUG_REFFRAME >=2
   cerr << "+++ void tle(ivg::Crf * crf_ptr,  const string path)" << endl; 
   tictoc tim;
   tim.tic();
#endif
   
    // assignment between PRN/PG and NORAD catalog number
    map<int,string> norad = {{28474,"PG02"}, {38833,"PG24"}, {29601,"PG12"}, {36585,"PG25"}};
    
    ifstream inStream;
    string line;

    while( ivg::parser::get_line(path,inStream, line))
    {
        if(line.substr(0,1) != "#")
        {
            string line1 = line;
            string line2;
            ivg::parser::get_line(path,inStream, line2);

            tle_t tle;
            int err_code = parse_elements( line1.c_str(), line2.c_str(), &tle);

            ivg::Source * source;
            if(crf_ptr->get_source( &source, norad[tle.norad_number] ))
                source->set_tle(tle);
        }
    }
    
    inStream.close();
    
    log<INFO>("*** ivg::parser::tle Loaded TLE-elements from ") % path;
   
#if DEBUG_REFFRAME >=2
    cerr << "--- void void tle(ivg::Crf * crf_ptr,  const string path)" << " : " << tim.toc() << " s " << endl;
#endif 
}
// ...........................................................................
istream& get_line(const string &path, ifstream  &inStream, string &line)
// ...........................................................................
{
    #define TRY_TO_OPEN 5
    
    if( inStream.is_open() )
        return(getline( inStream, line ));
    else
    {
        for(int i=0; i<=TRY_TO_OPEN; i++)
        {
            inStream.open(path.c_str(), ios::in);
            if(inStream.is_open())
            {
                return(getline( inStream, line ));
            }
            else if(!inStream.is_open() && i < TRY_TO_OPEN)
            {
                log<INFO>("*** Waiting 1sec, trying [") % i % "/" % TRY_TO_OPEN % "] to open file: " % path;
                sleep(5);
            }
            else if( !inStream.is_open() && i == TRY_TO_OPEN)
            {
                stringstream errormessage;
                errormessage <<
                             "void get_line(const string path, ifstream  &inStream) in parser.cpp: Failed to open file: "
                             << path << endl;
                throw runtime_error( errormessage.str() );
            }
        }
        
        
        return(getline( inStream, line ));
    }

}
// ...........................................................................
int init_options(string option_str)
// ...........................................................................
{
#if DEBUG_REFFRAME >=2
   cerr << "+++ int init_options(string option_str)" << endl; 
   tictoc tim;
   tim.tic();
#endif
   
    int option = 0x0;
    stringstream ss(option_str);
    string s;

    while (getline(ss, s, '|')) {
        s = remove_spaces(s);
        
        if( s == "EST2TRF")
            option |= EST2TRF;
        else if( s == "EST2CRF")
            option |= EST2CRF;
        else if( s == "EST2EOP")
            option |= EST2EOP;
        else if( s == "APR2TRF")
            option |= APR2TRF;
        else if( s == "APR2CRF")
            option |= APR2CRF;
        else if( s == "APR2EOP")
            option |= APR2EOP;
        else if( s == "TRF2APR")
            option |= TRF2APR;
        else if( s == "CRF2APR")
            option |= CRF2APR;
        else if( s == "EOP2APR")
            option |= EOP2APR;
        else if( s == "NUT2ZER")
            option |= NUT2ZER;
        else
            log<WARNING>("!!! Ignoring unknown snx-adjustment-option: ") % s;
    }
    
#if DEBUG_REFFRAME >=2
    cerr << "--- int init_options(string option_str)" << " : " << tim.toc() << " s " << endl;
#endif 
    return option;
}
// ...........................................................................
} // # namespace parser
} // # namespace ivg
