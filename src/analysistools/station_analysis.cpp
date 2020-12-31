#include "station_analysis.h"

namespace ivgat
{


// ===========================================================================
//              constructors and destructor
// ===========================================================================
Station_analysis::Station_analysis( ){ }
// ...........................................................................
Station_analysis::Station_analysis( std::vector<ivg::Trf> trf_list )
// ...........................................................................
{
   _trf_list = trf_list;
   _use_trf = vector<bool>(_trf_list.size(),true);
   
   std::vector<ivg::Trf>::iterator it_trf;
   std::vector<ivg::Analysis_station>::iterator it_sta;
   
   for( it_trf = _trf_list.begin(); it_trf < _trf_list.end(); it_trf++ )
   {           
        for( it_sta = it_trf->begin(); it_sta < it_trf->end(); it_sta++ )
        {   
            _counter[ it_sta->get_name(ivg::staname::ivs_name) ]++;
        }
   }  

}

// ...........................................................................
void Station_analysis::check_stations( vector<ivg::Trf> trf0, double th )
// ...........................................................................
{
#if DEBUG_ANALYSIS >= 2
   cerr << "+++ void Station_analysis::check_stations( vector<ivg::Trf>, double )" << endl; 
   tictoc tim;
   tim.tic();
#endif

   std::vector<ivg::Trf> trf_list_mod; 
   std::vector<ivg::Trf>::iterator it_trf, it_trf0;
   std::vector<ivg::Analysis_station>::iterator it_sta, it_sta0;
   bool is_trf = true;
                        
   for( it_trf = _trf_list.begin(); it_trf < _trf_list.end(); it_trf++ )
   {
       is_trf = true; 
       int idx = it_trf - _trf_list.begin();        
      
        for( it_trf0 = trf0.begin(); it_trf0 < trf0.end(); ++it_trf0 )
        {
            for( it_sta = it_trf->begin(); it_sta < it_trf->end(); it_sta++ )
            {         
                for( it_sta0 = it_trf0->begin(); it_sta0 < it_trf0->end(); it_sta0++ )
                {                
                    if( it_sta->get_name( ivg::staname::ivs_name ) == it_sta0->get_name( ivg::staname::ivs_name ) 
                        && it_sta->get_refepoch() == it_sta0->get_refepoch() )
                    {                                    
            
                        ivg::Matrix diff = it_sta->get_xyz0() - it_sta0->get_xyz0();
                        ivg::Matrix d = ( it_sta->form_topo2geo() ).transpose() * diff;

                        if( (d.absD())(0,0) >= th || (d.absD())(1,0) >= th || (d.absD())(2,0) >= th )
                        {
                           is_trf = false;
                        }                       

                    }         
                }      
            }   
        } 
        if( is_trf == true )
            trf_list_mod.push_back( *it_trf );             
   }    
   
   if( _trf_list.size()-trf_list_mod.size() != 0 )
   {
      log<RESULT>("*** ") % (_trf_list.size()-trf_list_mod.size()) % " sessions removed from solution. Too large estimates for station position!";
      _trf_list = trf_list_mod; 
   }
   
#if DEBUG_ANALYSIS >= 2
   cerr << "--- void Station_analysis::check_stations( vector<ivg::Trf>, double )"
        << " : " << tim.toc() << " s " << endl; 
#endif   
}


// ...........................................................................
void Station_analysis::get_station_position( std::string sta, std::string type, ivg::Matrix & sta_ts, 
                                             ivg::Matrix & sta_std_ts, ivg::Matrix & epo_ts )
// ...........................................................................
{
#if DEBUG_ANALYSIS >= 3
   cerr << "+++ void Station_analysis::get_station_positions()" << endl; 
   tictoc tim;
   tim.tic();
#endif
         
    sta_ts.resize(3,_stations[sta].size(),0.0);
    sta_std_ts.resize(3,_stations[sta].size(),0.0);
    epo_ts.resize(1,_stations[sta].size(),0.0);
       
    for( int i=0; i<_stations[sta].size(); ++i )
    {
        if( type == "REN")
        {
            sta_ts.set_sub( 0,i,_stations[sta].at(i).ren );
            // sta_std_ts.set_sub(0,i,_stations[sta].at(i).ren_std );
        }
        else if( type == "XYZ")
        {
            sta_ts.set_sub( 0,i,_stations[sta].at(i).xyz );
            sta_std_ts.set_sub( 0,i,_stations[sta].at(i).xyz_std );
        }
        else
        {
            throw runtime_error( "Unkwon coordinate type selected. Use either REN or XYZ.\n" );
        }   
        epo_ts(0,i) = _stations[sta].at(i).epoch.get_double_mjd();
    }
         
#if DEBUG_ANALYSIS >= 3
   cerr << "--- void Station_analysis::get_station_positions()"
        << " : " << tim.toc() << " s " << endl; 
#endif

   return;
}


// ...........................................................................
void Station_analysis::calc_station_positions( vector<ivg::Trf> trf0, double th )
// ...........................................................................
{
#if DEBUG_ANALYSIS >= 1
   cerr << "+++ void Station_analysis::get_station_positions()" << endl; 
   tictoc tim;
   tim.tic();
#endif
   
   // right now "_use_trf" is not really used. It's just the first step towards a time optimization
   // we still need "check_stations" beforce "calc_station_positions"
   
   station_epochs sta_epo;
   std::vector<ivg::Trf> trf_list_mod; 
   std::vector<ivg::Trf>::iterator it_trf, it_trf0;
   std::vector<ivg::Analysis_station>::iterator it_sta, it_sta0;
   
   for( it_trf = _trf_list.begin(); it_trf < _trf_list.end(); it_trf++ )
   {
        int idx = it_trf - _trf_list.begin();   
      
        for( it_trf0 = trf0.begin(); it_trf0 < trf0.end(); ++it_trf0 )
        {
            if(_use_trf.at(idx) == false)
                break;
            
            for( it_sta = it_trf->begin(); it_sta < it_trf->end(); it_sta++ )
            {
                if(_use_trf.at(idx) == false)
                    break;
                
                for( it_sta0 = it_trf0->begin(); it_sta0 < it_trf0->end(); it_sta0++ )
                {       
                    if( it_sta->get_name( ivg::staname::ivs_name ) == it_sta0->get_name( ivg::staname::ivs_name ) 
                        && it_sta->get_refepoch() == it_sta0->get_refepoch() )
                    {                

                        ivg::Date epoch = it_sta->get_refepoch();                                               
                        ivg::Matrix xyz = it_sta->get_xyz0() - it_sta0->get_xyz0();
                        ivg::Matrix xyz_std = it_sta->get_xyz0_std();
                        // xyz = it_sta->get_xyz0() - ( it_sta0->get_xyz0() + it_sta0->calc_psd_displacement( epoch ) );
                        
                        ivg::Matrix ren = ( it_sta->form_topo2geo() ).transpose() * xyz;                        
                        // to-do: VFG for REN 
                        
                        sta_epo = { ren, xyz, xyz_std, epoch, idx };
                        
                        if( (ren.absD())(0,0) < th && (ren.absD())(1,0) < th && (ren.absD())(2,0) < th )
                            _stations[ it_sta->get_name( ivg::staname::ivs_name ) ].push_back( sta_epo ); 
                        else
                            _use_trf.at(idx) = false;
                              
                        break;                                        
                    }         
                }      
            }       
        }
       
        if( _use_trf.at(idx) == true )
            trf_list_mod.push_back( *it_trf );    
   }
   
    if( _trf_list.size()-trf_list_mod.size() != 0 )
    {
        log<RESULT>("*** ") % (_trf_list.size()-trf_list_mod.size()) % " sessions removed from solution. Too large estimates for station position!";
        _trf_list = trf_list_mod; 
    }

#if DEBUG_ANALYSIS >= 1
   cerr << "--- void Station_analysis::get_station_positions()"
        << " : " << tim.toc() << " s " << endl; 
#endif

   return;
}


// ...........................................................................
void Station_analysis::get_level_of_uncertainties( std::string name, std::string type,
                                                   ivg::Matrix & xyz_std, ivg::Matrix & epo )
// ...........................................................................
{
#if DEBUG_ANALYSIS >= 1
   cerr << "+++ void Station_analysis::get_station_positions()" << endl; 
   tictoc tim;
   tim.tic();
#endif
      
   std::vector<ivg::Trf>::iterator it_trf;
   std::vector<ivg::Analysis_station>::iterator it_sta;
   xyz_std.resize( 3, _counter[name]+1, 0.0 );
   std::vector<double> mjd;
   int c = 0;
   
   for( it_trf = _trf_list.begin(); it_trf < _trf_list.end(); it_trf++ )
   {           
        ivg::Analysis_station * sta; 
        bool found = it_trf->get_station( &sta, name ); 
                
        if( found )
        {
            mjd.push_back( it_trf->get_reference_epoch().get_double_mjd() );
            
            xyz_std( 0,c ) = sta->get_xyz0_std()(0);
            xyz_std( 1,c ) = sta->get_xyz0_std()(1);
            xyz_std( 2,c ) = sta->get_xyz0_std()(2);
            c++;
        }    
   }
   ivg::Matrix tmp( mjd );   
   epo = tmp;

#if DEBUG_ANALYSIS >= 1
   cerr << "--- void Station_analysis::get_station_positions()"
        << " : " << tim.toc() << " s " << endl; 
#endif

   return;
}


// ...........................................................................
void Station_analysis::calc_baseline_length()
// ...........................................................................
{
#if DEBUG_ANALYSIS >= 1
   cerr << "+++ void Station_analysis::calc_baseline_length()" << endl; 
   tictoc tim;
   tim.tic();
#endif

   double bl_length, bl_length_std;
   ivg::Matrix bl(1,3,0.0);
   std::string bl_name;
   std::vector<ivg::Trf>::iterator it_trf;
   std::vector<ivg::Analysis_station>::iterator it_stai, it_staj;
   
   for( it_trf = _trf_list.begin(); it_trf < _trf_list.end(); it_trf++ )
   {
      std::map<std::string, ivg::Matrix > bl_info;
      
      for( it_stai = it_trf->begin(); it_stai < it_trf->end(); it_stai++ )
      {
         
         ivg::Date epoch = it_stai->get_refepoch(); 
          
         for( it_staj = it_stai+1; it_staj < it_trf->end(); it_staj++ )
         {
            // get baseline names (e.g., "WETTZELL-WESTFORD") and calculate baseline lengths and 
            // save them to a map (together with the corresponding epoch) 
            bl_name = it_stai->get_name( ivg::staname::ivs_name )+"-"+it_staj->get_name( ivg::staname::ivs_name );
            bl_length = it_stai->calc_baseline_length( *it_staj, epoch, bl_length_std ); 

            bl(0,0) = epoch.get_double_mjd(); 
            bl(0,1) = bl_length;
            bl(0,2) = bl_length_std;
            bl_info[ bl_name ] = bl ;
            
            if( std::find(_all_bls.begin(), _all_bls.end(), bl_name) == _all_bls.end() )
                _all_bls.push_back( bl_name );     
         }
      }

      _baselines.push_back( bl_info );
   }
   
   sort( _all_bls.begin(), _all_bls.end() );   

#if DEBUG_ANALYSIS >= 1
   cerr << "--- void Station_analysis::calc_baseline_length()"
        << " : " << tim.toc() << " s " << endl; 
#endif

   return;
}


// ...........................................................................
void Station_analysis::show_baselines()
// ...........................................................................
{
#if DEBUG_ANALYSIS >= 3
   cerr << "+++ void Station_analysis::show_baselines()" << endl; 
   tictoc tim;
   tim.tic();
#endif

   typedef std::map<std::string, ivg::Matrix>::const_iterator mapIter; 
   typedef std::vector< std::map<std::string, ivg::Matrix> >::const_iterator vecIter;

   cerr << endl;
   cerr << setiosflags(ios::left) << "Station_analysis::show_baselines() " << endl;
   cerr << " " << setfill(' ') << setw(20) << "BASELINE"
        << " " << setfill(' ' ) << setw(15) << "MJD" 
        << " " << setfill(' ' ) << setw(15) << "BASELINE LENGTHS" 
        << " " << setfill(' ' ) << setw(15) << "BASELINE LENGTHS STDDEV"  << endl; 

   for( vecIter iter = _baselines.begin(); iter != _baselines.end(); iter++ )
   {
      for (mapIter it = (*iter).begin(); it != (*iter).end(); it++)
      {
         cerr << setiosflags(ios::left) << " " << setfill(' ') << setw(20) << it->first
              << " " << setfill(' ' ) << setw(15) << setprecision(12) << it->second(0,0) 
              << " " << setfill(' ' ) << setw(15) << setprecision(12) << it->second(0,1)
              << " " << setfill(' ' ) << setw(15) << setprecision(12) << it->second(0,2)   
              << endl; 
      }
   }

#if DEBUG_ANALYSIS >= 3
   cerr << "--- void Station_analysis::show_baselines()" 
        << " : " << tim.toc() << " s " << endl; 
#endif

   return;
}


// ...........................................................................
void Station_analysis::calc_bl_rep( std::vector<std::string> & bl_rep_name, ivg::Matrix & bl_rep, int n_bl_length, std::string type )
// ...........................................................................
{
#if DEBUG_ANALYSIS >= 1
   cerr << "+++ void Station_analysis::calc_bl_rep()" << endl; 
   tictoc tim;
   tim.tic();
#endif

    int nobs = 2000;
    int counter = 0;
    std::vector<std::string> bl_done;

    ivgat::Fit fitter;
    ivg::Matrix bl_tmp(1,2,0.0);

    typedef std::vector< std::map<std::string, ivg::Matrix> >::iterator vecIter;
     
    for( std::vector< std::string >::iterator it = _all_bls.begin(); it != _all_bls.end(); it++ )
    {
        ivg::Matrix bl_length( nobs, 1, 0.0 );
        ivg::Matrix bl_length_std( nobs, 1, 0.0 );
        ivg::Matrix mjd( nobs, 1, 0.0 );
        std::vector<std::string> bl_name;
        counter = 0;

        for( vecIter iter2 = _baselines.begin(); iter2 != _baselines.end(); iter2++ )
        {
             auto search = (*iter2).find( *it );

             if( search != (*iter2).end() )
             {
                bl_done.push_back( search->first );
                bl_name.push_back( search->first );
                bl_length_std( counter, 0 ) = search->second(0,2);
                bl_length( counter, 0 ) = search->second(0,1);
                mjd( counter, 0 ) = search->second(0,0);
                counter++;
             }
        }
      
        std::vector<int> find_null = bl_length.find_idx( 0.0 );
        bl_length.rem_r( find_null );
        bl_length_std.rem_r( find_null );
        mjd.rem_r( find_null );

//        ivg::Matrix out = mjd;
//        out.append_cols(bl_length);
//        
//        std::string bas = "/home/halsig/ascot/output/bl_length/"+(*it)+".bin";
//        out.save_bin(bas);        
        
        if( bl_length.rows() >= n_bl_length )
        {
           ivg::Matrix bs = fitter.polyfit( mjd, bl_length, 1.0 );
           
           // get post-fit residuals
           ivg::Matrix r = fitter.get_resid();

           double rms;
           if( type == "RMS" )
               rms = sqrt( ( r.transpose()*r )(0) / r.size(1) );
           else if( type == "WRMS" )
           {                       
                ivg::Matrix p;
                p.eye( bl_length_std.rows() );
                p.solve_neq( bl_length_std.diag() );

                double sum = 0.0;
                double sum_p = 0.0;
                double tmp;
                for( int k=0; k<r.rows(); k++ )
                {
                    tmp = 1.0 / pow( bl_length_std(k), 2.0 );
                    sum_p += tmp;
                    sum += pow( r(k),2.0 ) * tmp;
                }
                rms = sqrt(sum / sum_p); 
           }

           bl_tmp(0,0) = bs.meanD();
           bl_tmp(0,1) = rms;
           _bl_rep[ *it] = bl_tmp;       
        }
    }   
      
#if DEBUG_ANALYSIS >= 1
   cerr << "--- void Station_analysis::calc_bl_rep()"
        << " : " << tim.toc() << " s " << endl; 
#endif

   return;
}


// ...........................................................................
void Station_analysis::show_bl_repeatabilities()
// ...........................................................................
{
#if DEBUG_ANALYSIS >= 3
   cerr << "+++ void Station_analysis::show_bl_repeatabilities()" << endl; 
   tictoc tim;
   tim.tic();
#endif
   
   typedef std::map<std::string, ivg::Matrix>::const_iterator mapIter; 

   cerr << endl;
   cerr << setiosflags(ios::left) << "Station_analysis::show_bl_repeatabilities() " << endl;
   cerr << " " << setfill(' ') << setw(20) << "BASELINE NAME"
        << " " << setfill(' ' ) << setw(17) << "BASELINE LENGTHS" 
        << " " << setfill(' ' ) << setw(17) << "B.L. REPEATABILITIES" << endl;
   cerr << "=======================================================" << endl; 

   for (mapIter it = _bl_rep.begin(); it != _bl_rep.end(); it++)
   {
      cerr << setiosflags(ios::left) << " " << setfill(' ') << setw(20) << it->first
           << " " << setfill(' ' ) << setw(17) << setprecision(12) << it->second(0,0) 
           << " " << setfill(' ' ) << setw(17) << setprecision(8) << it->second(0,1)
           << endl; 
   }

   cerr << "=======================================================" << endl;
   cerr << "#baselines: " << _bl_rep.size() << endl;
   
#if DEBUG_ANALYSIS >= 3
   cerr << "--- void Station_analysis::show_bl_repeatabilities()" 
        << " : " << tim.toc() << " s " << endl; 
#endif

   return;
}


// ...........................................................................
void Station_analysis::remove_baselines( std::vector<std::string> rm_bls )
// ...........................................................................
{
#if DEBUG_ANALYSIS >= 1
   cerr << "+++ void remove_baselines( std::vector<std::string> )" << endl; 
   tictoc tim;
   tim.tic();
#endif
   
    std::map<std::string, ivg::Matrix>::iterator iter;

    for( std::vector< std::string >::iterator it = rm_bls.begin(); it != rm_bls.end(); it++ )
    {
        iter = _bl_rep.find( *it );
       
        if( iter != _bl_rep.end() )
            _bl_rep.erase( iter );
    }

#if DEBUG_ANALYSIS >= 1
   cerr << "--- void remove_baselines( std::vector<std::string> )" 
        << " : " << tim.toc() << " s " << endl; 
#endif

   return;
}


// ...........................................................................
std::map<std::string, ivg::Matrix> Station_analysis::fit_bl_rep( ivg::Matrix & bl )
// ...........................................................................
{
#if DEBUG_ANALYSIS >= 2
   cerr << "+++ void Station_analysis::fit_bl_rep()" << endl; 
   tictoc tim;
   tim.tic();
#endif

   bl.resize( _bl_rep.size(),2,0.0 );
   int c = 0;
   std::map<std::string, ivg::Matrix> fit;
   fit = _bl_rep;
   typedef std::map<std::string, ivg::Matrix>::const_iterator mapIter; 

   for (mapIter it = _bl_rep.begin(); it != _bl_rep.end(); it++)
   {
      bl( c,0 ) = it->second(0,0);
      bl( c,1 ) = it->second(0,1);
      c++;
   }

   // fit trend (e.g., quadratic polynomial) to the baseline lengths 
   // repeatabilities and evaluate for each baseline lengths
   ivgat::Fit fitter;
   //ivg::Matrix bs = fitter.polyfit( bl.get_col(0), bl.get_col(1), 2.0 );
   ivg::Matrix bs = fitter.expfit( bl.get_col(0), bl.get_col(1) );
   bl.set_sub(0,1,bs);
   
   c = 0;
   for (mapIter it = _bl_rep.begin(); it != _bl_rep.end(); it++)
   {
      fit[ it->first ] = bl.get_row(c);
      c++;
   }
  
#if DEBUG_ANALYSIS >= 2
   cerr << "--- void Station_analysis::fit_bl_rep()"
        << " : " << tim.toc() << " s " << endl; 
#endif

   return fit;
}

// ...........................................................................
void Station_analysis::calc_RMS_differences( std::map< std::string, ivg::Matrix > other, 
                                             std::string sol1, std::string sol2,
                                             double amount )
// ...........................................................................
{
#if DEBUG_ANALYSIS >= 3
   cerr << "+++ void calc_RMS_differences( std::map< std::string, ivg::Matrix > )" << endl; 
   tictoc tim;
   tim.tic();
#endif

    ivg::Matrix d_bl_rep( _bl_rep.size(),1,0.0 );
    int c = 0;

    typedef std::map<std::string, ivg::Matrix>::const_iterator mapIter; 
    mapIter it_other;

    for (mapIter it = _bl_rep.begin(); it != _bl_rep.end(); it++)
    {
        it_other = other.find( it->first );

        if( it_other != other.end() )
        {          
             d_bl_rep(c) = it->second(0,1) - it_other->second(0,1);
             c++;
        }         
    }   

    std::vector<int> idx1 = d_bl_rep.find_idx( lt, -amount );
    std::vector<int> idx2 = d_bl_rep.find_idx( gt, amount );
    
    double improve = 100.0 / d_bl_rep.rows() * idx1.size();
    double degrade = 100.0 / d_bl_rep.rows() * idx2.size();
    double unchanged = 100.0 - improve - degrade ;

    log<NOTHING>("*** Solution '") % sol2 % "' w.r.t. Solution '" % sol1 
        % "': improvement of at least 1mm for " % improve 
        % "[%], degradation of at least 1mm for: " % degrade  
        % ": [%] of baselines; " % unchanged % "[%] of baselines remain unchanged.";
   
   
#if DEBUG_ANALYSIS >= 3
   cerr << "--- void calc_RMS_differences( std::map< std::string, ivg::Matrix > )"
        << " : " << tim.toc() << " s " << endl; 
#endif

}

// ...........................................................................
ivg::Matrix Station_analysis::calc_helmert_params( vector<ivg::Trf> other, vector<t_param> tp )
// ...........................................................................
{
#if DEBUG_ANALYSIS >= 3
   cerr << "+++ void Station_analysis::calc_helmert_params()" << endl; 
   tictoc tim;
   tim.tic();
#endif
      
   std::vector<ivg::Trf>::iterator it_trf;
   std::vector<ivg::Analysis_station>::iterator it_sta;
   ivg::Matrix result_series;
   ivgat::Transformation trans( ivgat::t_type::cart_3D );
   
   for( it_trf = _trf_list.begin(); it_trf < _trf_list.end(); it_trf++ )
   {
       
        ivg::Date epoch = it_trf->get_station(0)->get_refepoch();
        ivg::Matrix trf( it_trf->get_number_stations(), 3, 0.0 ); 
        ivg::Matrix trf2( it_trf->get_number_stations(), 3, 0.0 ); 

        ivg::Trf tmp = other.at( it_trf-_trf_list.begin() );
        ivg::Matrix idx_sta = it_trf->get_corresponding_stations( tmp );
        int k = 0; 
        
        for( it_sta = it_trf->begin(); it_sta < it_trf->end(); it_sta++ )
        {
            trf.set_sub( k,0, it_sta->get_xyz0().transpose() );
            k++;
        }  

        for( int i = 0; i < idx_sta.size(1); i++ )
        {
            trf2.set_sub( i,0, tmp.get_station(i)->get_xyz0().transpose() );
        }

        trans.set_systems( trf, trf2 );

        // define return matrices (call by reference) and the parameter which should be estimated
        ivg::Matrix hparams, std, corr;      
        if( trf.size(1) > tp.size() )
        {
            trans.estimate_parameter( tp, hparams, std, corr );
            ivg::Matrix tmp(1,1,epoch.get_double_mjd());
            tmp.append_rows(hparams);
            // containing [mjd, #params]
            result_series.append_rows(tmp.transpose());
        }  
   }

 
#if DEBUG_ANALYSIS >= 3
   cerr << "--- void Station_analysis::calc_helmert_params()"
        << " : " << tim.toc() << " s " << endl; 
#endif

   return result_series;
}



// ...........................................................................
std::vector<std::string> Station_analysis::get_sessions()
// ...........................................................................
{
#if DEBUG_ANALYSIS >= 3
   cerr << "+++ std::vector<std::string> ivg::Station_analysis get_sessions()" << endl; 
   tictoc tim;
   tim.tic();
#endif

   std::vector<ivg::Trf>::iterator it_trf;
   std::vector<std::string> sessions;
   
   for( it_trf = _trf_list.begin(); it_trf < _trf_list.end(); it_trf++ )
   {      
        std::string s = it_trf->get_name();
        std::string delim = "/";

        auto start = 0U;
        auto end = s.find(delim);
        while (end != std::string::npos)
        {
            start = end + delim.length();
            end = s.find(delim, start);
        }
       
        sessions.push_back( s.substr(start, end) );
   }   

#if DEBUG_ANALYSIS >= 3
   cerr << "--- std::vector<std::string> ivg::Station_analysis get_sessions()"
        << " : " << tim.toc() << " s " << endl; 
#endif
   
   return sessions;
}


} // # namespace ivgat
