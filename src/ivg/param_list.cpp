#include "param_list.h"

namespace ivg
{

// ...........................................................................
Param_list::Param_list()
// ...........................................................................
{
}

// ...........................................................................
Param_list::~Param_list()
// ...........................................................................
{
}

// ...........................................................................
Param_list::Param_list(const Param_list& orig)
// ...........................................................................
{
    // parameter vector
    _params.resize( orig._params.size() );
    copy( orig._params.begin(), orig._params.end(),_params.begin() );
    
    _start_epoch = orig._start_epoch;
    _end_epoch = orig._end_epoch;
}

// ...........................................................................
Param_list::Param_list( const vector<ivg::Param> params )
// ...........................................................................
{
   _params = params;
}

// ...........................................................................
Param_list::Param_list( const vector<ivg::Param> params, const vector<ivg::Param> stoch_params )
// ...........................................................................
{
   _params = params;
   _stoch_params = stoch_params;
   
   double min_mjd = 10.0e10;
   double max_mjd = 0.0;
   for( auto & param : _stoch_params )
   {
       double mjd = param.get_epoch().get_double_mjd();
       
       if( mjd < min_mjd )
           min_mjd = mjd;
       if( mjd > max_mjd )
           max_mjd = mjd;
   }
   _start_epoch = ivg::Date( min_mjd );
   _end_epoch = ivg::Date( max_mjd );
}

// ...........................................................................
Param_list::Param_list(ivg::Trf &trf, ivg::Crf &crf, ivg::Eop_series &eops, ivg::Date epoch)
// ...........................................................................
{

    for( vector<ivg::Analysis_station>::iterator sta_iter = trf.begin();
            sta_iter!=trf.end(); ++sta_iter )
    {
        string sta_name = sta_iter->get_name( ivg::staname::ivs_name );
        // calculate station position at reference epoch; use PSD in every case
        // as it is only initialized if it is also applied
        ivg::Matrix sta_xyz = sta_iter->calc_xyz(epoch,{"PSD"})/ivg::c*1e2; // m -> centi-seconds

        _params.push_back(ivg::Param(ivg::paramtype::stax, sta_name, epoch, sta_xyz(0)));
        _params.push_back(ivg::Param(ivg::paramtype::stay, sta_name, epoch, sta_xyz(1)));
        _params.push_back(ivg::Param(ivg::paramtype::staz, sta_name, epoch, sta_xyz(2)));
        _params.push_back(ivg::Param(ivg::paramtype::clo, sta_name, epoch));
        _params.push_back(ivg::Param(ivg::paramtype::zwd, sta_name, epoch));
        _params.push_back(ivg::Param(ivg::paramtype::ngr, sta_name, epoch));
        _params.push_back(ivg::Param(ivg::paramtype::egr, sta_name, epoch));
    }

    ivg::Matrix erp = eops.calc_erp(epoch);
    ivg::Matrix nut = eops.calc_nut(epoch);
    _params.push_back(ivg::Param(ivg::paramtype::xpo, "EOP", epoch,erp(0)));
    _params.push_back(ivg::Param(ivg::paramtype::ypo, "EOP", epoch,erp(1)));
    _params.push_back(ivg::Param(ivg::paramtype::ut1, "EOP", epoch,erp(2)));
    _params.push_back(ivg::Param(ivg::paramtype::nutx, "EOP", epoch,nut(0)));
    _params.push_back(ivg::Param(ivg::paramtype::nuty, "EOP", epoch,nut(1)));


    for( vector<ivg::Source>::iterator src_iter = crf.begin();
            src_iter!=crf.end(); ++src_iter )
    {
        if(src_iter->use_me())
        {
            string src_name = src_iter->get_name( ivg::srcname::iers );
            _params.push_back(ivg::Param(ivg::paramtype::ra, src_name, epoch,src_iter->get_ra0()));
            _params.push_back(ivg::Param(ivg::paramtype::dec, src_name, epoch,src_iter->get_dec0()));
        }
    }
     
    // baseline-dependent clock offsets (baseline clocks)
    vector<string> lc = trf.get_station_names(ivg::staname::lettercode);
    for( int i=0; i<lc.size(); ++i )
    {       
        for( int j=i+1; j<lc.size(); ++j )
        {
            string bl = lc.at(i)+"-"+lc.at(j);
            _params.push_back(ivg::Param(ivg::paramtype::blcl,bl,epoch));
        }
    }    
}


// ...........................................................................
void Param_list::modify_parameterization( std::string dbname, Setting &setup, 
                                          ivg::Trf &trf, ivg::Crf &crf, 
                                          ivg::Ls_solution &solver, ivg::Matrix apriori,
                                          ivg::Matrix obs_epochs )
// ...........................................................................
{
#if DEBUG_VLBI >=2
   cerr << "+++ void Param_list::modify_parameterization(string, Setting &, "
        << "ivg::Trf &, ivg::Crf &, ivg::Ls_solution & )" << endl; 
   tictoc tim;
   tim.tic();
#endif
    
    ivg::Ls_solution* sol = &solver;
    ivg::Ls_neq* neq = dynamic_cast<ivg::Ls_neq*>(sol);
    ivg::Lsa* gmm    = dynamic_cast<ivg::Lsa*>(sol);
    ivg::Lsc* lsc    = dynamic_cast<ivg::Lsc*>(sol);
       
    // insert clock breaks
    if( gmm ) 
    {
        _insert_breaks( gmm, ivg::paramtype::cbr );
        _insert_breaks( gmm, ivg::paramtype::atbr );
    }
    else if( lsc )
    {
        _insert_breaks( lsc->get_lsa_ptr(), ivg::paramtype::cbr );
        _insert_breaks( lsc->get_lsa_ptr(), ivg::paramtype::atbr );        
    }
   
    //Extracting PARAMS group from config file
    Setting &params = setup[ "PARAMS" ];

    map<string, vector<ivg::paramtype> > assignment =
    {
        { "stations", {ivg::paramtype::stax, ivg::paramtype::stay, ivg::paramtype::staz}},
        { "clocks", {ivg::paramtype::clo,ivg::paramtype::cbr}} ,
        { "zwd", {ivg::paramtype::zwd, ivg::paramtype::atbr}} ,
        { "ut1", {ivg::paramtype::ut1}} ,
        { "pm", {ivg::paramtype::xpo,ivg::paramtype::ypo}} ,
        { "nut", {ivg::paramtype::nutx,ivg::paramtype::nuty}} ,
        { "gradients", {ivg::paramtype::ngr,ivg::paramtype::egr}} ,
        { "sources", {ivg::paramtype::ra,ivg::paramtype::dec}},
        { "bl_clocks", {ivg::paramtype::blcl}}
    };

    // only if nutobln exists in the config-file, take it into account (for sinex combination dgfi)
    if(params.exists("nutobln"))
        assignment.insert({ "nutobln", {ivg::paramtype::nutob,ivg::paramtype::nutln}});

    // create list of        ...
    vector<string>sta_names = trf.get_station_names( ivg::staname::ivs_name );
    vector<string>src_names = crf.get_source_names( ivg::srcname::iers );
    vector<string>names;
    // baseline-dependent clock offsets (baseline clocks)
    vector<string> bl_names;
    vector<string> bl_altnames;
    vector<string> lc = trf.get_station_names(ivg::staname::lettercode);
    for( int i=0; i<lc.size(); ++i )
      for( int j=i+1; j<lc.size(); ++j ) {
            bl_names.push_back( lc.at(i)+"-"+lc.at(j) );
	    bl_altnames.push_back( sta_names.at(i)+"-"+sta_names.at(j) );
            bl_altnames.push_back( sta_names.at(j)+"-"+sta_names.at(i) );
	    }
    
    // loop over all necessary config files entries in PARAM-block
    for( map<string, vector<ivg::paramtype> >::iterator param_iter = assignment.begin(); param_iter != assignment.end(); ++param_iter)
    {
         vector<string> altnames;
        // select relevant param names
        if( param_iter->first == "stations" || param_iter->first == "zwd" ||
                param_iter->first == "gradients"  || param_iter->first == "clocks" )
	  {
            names = sta_names;
	   
	  }
        else if(param_iter->first == "sources")
	  {
	   names = src_names;
	   
	  }
        else if(param_iter->first == "bl_clocks")
	  {
            names = bl_names;
	    altnames = bl_altnames;
	  }
        else if(param_iter->first == "nut" || param_iter->first == "nutobln" || param_iter->first == "pm"
                || param_iter->first == "ut1" )
        {
            names.clear();
            names.push_back( "EOP" );
	    
        }
	vector<string> vgosdb_lst;
        // loop over list of this config-file entry
       
	//for( int j=0; j < params[ param_iter->first ].getLength(); ++j )
	for( int j=0; j < params.lookup( param_iter->first ).getLength(); ++j )
        {
	    
            string selected_group_name = params.lookup( param_iter->first )[ j ][0].getName(); // here  e.g. "stations"
            string handling = params.lookup( param_iter->first )[ j ]["handling"];

           // detect if "stations" is a string (e.g. vgosdb) then getType == 4 (stands for string)
           int name_idx;
           if( params.lookup( param_iter->first )[ j ][0].getType() == 4) // j == 0 && param_iter->first == "clocks" && 
           {
               string string_group = params.lookup( param_iter->first )[ j ][0]; // 0 is the item "stations"
               if( string_group == "vgosdb")
                    name_idx = -999;
               else
                   throw runtime_error( "!!! Error in config-file: In "+param_iter->first+", stations = only can be a number or exactly the string vgosdb. Nothing else.\n" );
           }
           else
                name_idx = params.lookup( param_iter->first )[ j ].lookup( selected_group_name );
            
            // select all param names of this class (stations, sources, ... ) within this session
            // group index 0  -> all stations, source or what ever
            //             >0 -> this group
            //             <0 -> all but this group
            //            -999 -> in case of group from vgosDB RefClockStationList 
            vector<string> param_lst;
            if( name_idx == -999)
            {
	      if(param_iter->first == "bl_clocks")
		{
		  if(setup.exists("BaselineClockList"))
                    {
		      for( int i=0; i<setup["BaselineClockList"].getLength(); ++i )
			{
     			  param_lst.push_back( setup["BaselineClockList"][ i ] );
			}
		    }
		  for( int i=0;i<_bl_clocks.size();i++)
		    param_lst.push_back( _bl_clocks.at(i));
		} else
		{
		  if(!setup.exists("RefClockStationList"))
                    throw logic_error( "!!! Wrong group definition (vgosdb) for reference clock in configfile. Reference Clock from vgosDB not set because processing non vgosdb based session.\n" );
		  
		  param_lst = {setup["RefClockStationList"]};
		}
	      vgosdb_lst=param_lst;
	      
            }else if(name_idx == 0)
            {
                param_lst = names;
		for( int i=0; i<vgosdb_lst.size(); ++i )
                {
                    vector<string>::iterator iter = find( param_lst.begin(), param_lst.end(),
                                                          vgosdb_lst[ i ] );
                    if (iter != param_lst.end())
                        param_lst.erase( iter );
                }
            }
            // select only those param-names within the requested group
            else if( name_idx > 0 )
            {
	      for( int i=0; i<setup["groups"].lookup(selected_group_name)[ name_idx-1 ].getLength(); ++i ) {
                    param_lst.push_back( setup["groups"].lookup(selected_group_name)[ name_idx -1 ][ i ] );

	      }
            }
            // remove entries of given group from the entire list
            else if( name_idx < 0)
            {
                name_idx *= -1;
                param_lst = names;
                for( int i=0; i<setup["groups"].lookup(selected_group_name)[ name_idx-1].getLength(); ++i )
                {
                    vector<string>::iterator iter = find( param_lst.begin(), param_lst.end(),
                                                          (const char *)setup["groups"].lookup(selected_group_name)[ name_idx-1][ i ] );
                    if (iter != param_lst.end())
                        param_lst.erase( iter );
                }
            }

	    if (altnames.size()>0)
	      {
		for (int i=0;i<param_lst.size();i++)
		  {
		    for (int k=0;k<altnames.size();k+=2)
		       {
			 if ((param_lst.at(i)==altnames.at(k))||(param_lst.at(i)==altnames.at(k+1)))
			   {
			     param_lst[i]=names[k/2];
			   }
			   
		       }
		  }
		for (int i=0;i<vgosdb_lst.size();i++)
		  {
		    for (int k=0;k<altnames.size();k+=2)
		       {
			 if ((vgosdb_lst.at(i)==altnames.at(k))||(vgosdb_lst.at(i)==altnames.at(k+1)))
			   {
			     vgosdb_lst[i]=names[k/2];
			   }
			   
		       }
		  }
		
	      }
            // create list of indices which are affected by cnf-entry
            vector<int> cur_idx;
            for(int i=0; i<param_lst.size(); i++)
            {
	      
                for(auto &tmp_type: assignment[ param_iter->first ])
                {
                    // changed from get_indexes to ge_idx in order to get ALL indexes and not only the first one
                    vector<int> fixis = get_idx(tmp_type, param_lst.at(i) );
                    
                    // if fixis is empty, we need to store a -1 for not-found
                    if(fixis.empty())
                        fixis.push_back(-1);
                    
                    cur_idx.insert( cur_idx.end(), fixis.begin(), fixis.end() );
                }
                
                // in case of clocks, find positions of all clock breaks and add them to cur_idx
                if( param_iter->first == "clocks" )
                {
                   vector<int> idx = get_idx( ivg::paramtype::cbr,param_lst.at(i) );
                   if( idx.size() != 0 )
                      cur_idx.insert( cur_idx.end(), idx.begin(), idx.end() );
                }

                // sort index vector and remove duplicates as well as -1
                sort( cur_idx.begin(), cur_idx.end() );
                cur_idx.erase( unique( cur_idx.begin(), cur_idx.end() ), cur_idx.end() );
                if( cur_idx.at( 0 ) == -1 )
                    cur_idx.erase( cur_idx.begin() );
            }

            // reduce params (set flag which indicates that parameter should be reduced;
            // this will be done in Session::solve)
            if( handling.find("reduce") != string::npos )
            {
	      //   vector<int> org_cur_idx=cur_idx;
                // check if more orders need to be reduced and then fix all
                int order = (int)params.lookup( param_iter->first )[ j ]["polynom"][ "order" ];
		// if( order != 0 )
                //{
		//  cur_idx.clear();

		    //  vector<int> orders;
		//  for(int tmp_order=1; tmp_order<=order; tmp_order++) {
	          //     orders.push_back(tmp_order);

		      //for(auto &tmp_param : _params)
		//     for (int i=0;i<_params.size();i++)
		//	{
			  //   if(tmp_param.is_type(param_iter->second, orders))
			  //if ((_params.at(i).get_type()==param_iter->second)&&(_params.at(i).get_order()==tmp_order))
		//	  for (int j=0;j<(param_iter->second).size();j++){
		//	    if (_params.at(i).is_type((param_iter->second)[j],order))
		//	      {
		//		cur_idx.push_back(i);
		//	      }
		//	  }
		//	}
                //    }
		//}
		
		if( order != 0 )
                {
		  //for (int i=0;i<cur_idx.size();i++)
		  //  std::cout << cur_idx[i]<< " ";
		  //std::cout << endl;
                    cur_idx.clear();
   
                    vector<int> orders;
                    for(int tmp_order=0; tmp_order<=order; tmp_order++)
                        orders.push_back(tmp_order);

                    //for(auto &tmp_param : _params)
		    for (int i=0;i<_params.size();i++)  
                    {
		      // if(tmp_param.is_type(param_iter->second, orders))
		      if ((_params.at(i).is_type(param_iter->second, orders))&&(_params.at(i).is_type_name(param_iter->second, param_lst)))
                        {
			  //  cur_idx.push_back(get_index(tmp_param));
			   cur_idx.push_back(i);
                        }
		      
                    }
                }
		//for (int i=0;i<cur_idx.size();i++)
		//   std::cout << cur_idx[i]<< " ";
		//  std::cout << endl;
                for( int i=0; i<cur_idx.size(); ++i ) {
                    _params.at( cur_idx.at( i ) ).set_reduce_flag( true );
		   
		}
            }
                        
            // include stochastic parameters
            if( handling.find("stochastic") != string::npos)
            {
                // check if more orders need to be fixed and then fix all
                int order = (int)params.lookup( param_iter->first )[ j ]["polynom"][ "order" ];
                if( order != 0 )
                {
                    cur_idx.clear();

                    vector<int> orders;
                    for(int tmp_order=0; tmp_order<=order; tmp_order++)
                        orders.push_back(tmp_order);

                    for(auto &tmp_param : _params)
                    {
                        if(tmp_param.is_type(param_iter->second, orders))
                        {
                            cur_idx.push_back(get_index(tmp_param));
                        }
                    }
                }
                
                solver.handle_stoch_param( cur_idx );    

                // include stochastic parameters (one parameter for each observation)
                // into the parameter list
                for( int i=0; i<cur_idx.size();++i)
                {
                    for( int n=1; n <= lsc->get_n_stoch_params().at(i); ++n )
                    {
                        _stoch_params.push_back( _params.at( cur_idx.at( i ) ) );
                        _stoch_params.back().set_order( -1 ); 
                        ivg::Matrix sta_epochs = _stoch_params_epochs[ _stoch_params.back().get_name() ];
                        ivg::Date obs_epo( sta_epochs(n-1,0) );
                        _stoch_params.back().set_epoch( obs_epo );                    
                    }                   
                }
//                for( int i=cur_idx.size()-1; i>=0; --i )
//                {                
//                    // delete deterministic offset parameter from param_list and apriori matrix
//                    _params.erase( _params.begin()+cur_idx.at( i ) );
//                    apriori.rem_c( cur_idx.at(i) );
//                }
            }
            
            // fix params
            if( handling.find("fix") != string::npos  )
            {
                // check if more orders need to be fixed and then fix all
	        
	        int order = (int)params.lookup( param_iter->first )[ j ]["polynom"][ "order" ];
                if( order != 0 )
                {
                    cur_idx.clear();

                    vector<int> orders;
                    for(int tmp_order=0; tmp_order<=order; tmp_order++)
                        orders.push_back(tmp_order);

                    for(auto &tmp_param : _params)
                    {
                        if(tmp_param.is_type(param_iter->second, orders))
                        {
                            cur_idx.push_back(get_index(tmp_param));
                        }
                    }
                }
		
                solver.fix_param( cur_idx ); 

                apriori.rem_c( cur_idx );
                for( int i=cur_idx.size()-1; i>=0; --i )
                {
                    vector<ivg::Param>::iterator it_param = _params.begin()+cur_idx.at( i );
                    log<DETAIL>("*** Fixing parameter name: ") % it_param->get_name() % " type: " % it_param->get_typename() % " order: " % it_param->get_order();

                    // for clocks: check if there is a clock break
                    if( param_iter->first == "clocks" )
                    {
                       vector<int> idx = get_idx( ivg::paramtype::cbr,it_param->get_name() );
                       if( idx.size() != 0 )
                       {   
                          // there is a break => do not fix, but stop with a logic error
                          throw logic_error( "void Param_list::modify_parameterization( ... ): Station "+it_param->get_name()+" has a clock break and must not be fixed." );
                       }
                    }

                    _params.erase( it_param );
                }
		
            }
            else
            {
                // if handling does NOT include the keyword "keep" (e.g. "keep_fix")
                // polynoms and/or cpwlf are set up
                if(handling.find("keep") == string::npos)
                {
		    
                    // transform constant parameters to a polynomial and/or piece-wise linear
                    // representation
                    for( int i=0; i<cur_idx.size(); ++i )
                    {
		      //if(param_iter->first == "bl_clocks")
                        // set offset constraint of constant parameter
                        _params.at( cur_idx.at( i ) ).set_offset_cnstr_sigma( (
                                    double)params.lookup( param_iter->first )[ j ]["polynom"][ "cnstr" ][ 0 ] );

                        // insert polynomials
                        if( (int)params.lookup( param_iter->first )[ j ]["polynom"][ "order" ] != 0 )
                        {
                            int order = (int)params.lookup( param_iter->first )[ j ]["polynom"][ "order" ];

                            ivg::Matrix new_apriori(order+1,1,0.0);
                            
                            //solver.trafo_params2polynomial( order, cur_idx.at( i ), _start_epoch.get_double_mjd() );
                            ivg::Date epo = _params.at( cur_idx.at( i ) ).get_epoch();
			    
			    solver.trafo_params2polynomial( order, cur_idx.at( i ),
                                                            epo.get_double_mjd(), apriori,
                                                            new_apriori );
			   
                            _params.at(cur_idx.at(i)).set_apriori( new_apriori(0) );
			  
                            for( int d=1; d<=order; ++d )
                            {
                                _params.push_back( _params.at( cur_idx.at( i ) ) );
                                _params.back().set_order( d );
                                _params.back().set_apriori( new_apriori(d,0) );
                                _params.back().set_offset_cnstr_sigma( (double)
                                                                       params.lookup( param_iter->first )[ j ]["polynom"][ "cnstr" ][ d ] );
				
                            }
			    
                        }
                        // if only an offset should be parametrized, the apirori
                        // value for this offset should be set to the parameter list
                        // (in case of cpwlf these "original" offsets are eliminated )
                        else
                        {
                            ivg::Date epo = _params.at( cur_idx.at( i ) ).get_epoch();          
                            ivg::Matrix apr_col = apriori( ":",cur_idx.at(i) );
                            std::vector<int> ind = apr_col.find_idx( ne, 0.0 );
                            
                            ivg::Matrix new_apriori, new_epo;
                            if( ind.begin() != ind.end() )
                            {
                                new_apriori = apr_col.get_sub( ind, {0} );                            
                                new_epo = obs_epochs.get_sub( ind, {0} );
                            
                                ivg::Matrix ones(new_epo.rows(),1,1.0);
                                ivg::Matrix diff = ( new_epo - ones*epo.get_double_mjd() ).absD();

                                int obs_idx;
                                diff.min(obs_idx);

                                _params.at(cur_idx.at(i)).set_apriori( new_apriori( obs_idx,0 ) );
                            }
                        }
                    }
                    
                    //loop over all parameters (indices) which are affected by controlfile-entry
                    for( int i=cur_idx.size()-1; i>=0; --i )
                    {
		     
                        // insert CPWLFs
                        if( (bool)params.lookup( param_iter->first )[ j ]["cpwlf"][ "insert" ] )
                        {
                            // special consideration of breaks
                            // when clock breaks or CPWLF breaks for ZWDs have been inserted, they should not be transformed to 
                            // CPWLF, however, they have to be combined with clock offsets 

                            // But, this is not necessary if we are dealing with a NEQ-solution
                            if( !neq && ( param_iter->first == "zwd" || param_iter->first == "clocks" ) )
                            {
                                if(  _params.at( cur_idx.at( i ) ).get_type() == ivg::paramtype::cbr
                                     || _params.at( cur_idx.at( i ) ).get_type() == ivg::paramtype::atbr )
                                    continue;
                                else
                                {
                                    // combine first offset paramater and following breaks to one single parameter
                                    // (but do not remove any of them). The first one will then be transformed to 
                                    // CPWLF, while the other ones (i.e., the breaks) remain in the equation
                                    // system. 

                                    // get_idx will find only the offset and no further higher degree polynomial parameters
                                    vector<int> idx_br;
                                    if( param_iter->first == "clocks" )
                                        idx_br = get_idx( ivg::paramtype::cbr, _params.at( cur_idx.at( i ) ).get_name() );
                                    else if( param_iter->first == "zwd" )
                                        idx_br = get_idx( ivg::paramtype::atbr, _params.at( cur_idx.at( i ) ).get_name() );
                                    
                                    if( idx_br.size() > 0 )
                                    {
                                        if( gmm )                                  
                                            gmm->combine_params( idx_br, false );
                                        else if( lsc )
                                            lsc->get_lsa_ptr()->combine_params( idx_br, false );  
                                    }
                                }
                            }
                            // end: special consideration of breaks
                            
                            // determine interval length
                            double dt = (double)params.lookup( param_iter->first )[ j ]["cpwlf"][ "int_length" ]/(60.0*24.0);
                    
                            // uses this constructor:  Matrix(double start, double schrittweite, double ende, int cols);
                            ivg::Matrix t_cpwlf( _start_epoch.get_double_mjd(), dt, _end_epoch.get_double_mjd()+dt/1.01, 1 );
                            
                            // only two CPWLF parameters and observations 
                            // outside interval => increase interval to end of
                            // session
                            if( t_cpwlf.rows() == 2 && _end_epoch.get_double_mjd() > t_cpwlf(1) )
                            {
                                t_cpwlf(1) = _end_epoch.get_double_mjd();
                                if(i == 1)
                                    log<WARNING>("!!! Extending interval length for ")  % param_iter->first % " from: " % dt % " to "  % ( t_cpwlf(1)-t_cpwlf(0) );
                            }
                            
                            // adding new interval for observations behind the last old interval
                            double timediff = _end_epoch.get_double_mjd()-t_cpwlf(t_cpwlf.length()-1);
                            if( timediff > 0.0 )
                                t_cpwlf.append_rows(t_cpwlf(t_cpwlf.length()-1)+dt);
			    
                            // transform equation system
                            solver.trafo_params2cpwlf( cur_idx.at( i ), t_cpwlf );
			   
                            // transform apriori to CPWLF parameters
                            ivg::Matrix apr_col = apriori( ":",cur_idx.at( i ) );
			    
                            // delete zero-rows in apriori-vector and remove corresponding epochs of observations
                            std::vector<int> ind = apr_col.find_idx( ne, 0.0 );                            
			     
                            ivg::Matrix new_apriori( t_cpwlf.rows(),1,0.0 );
			    
		            
                            if( ind.begin() != ind.end() )
                            {
                                ivg::Matrix apr_col_nz = apr_col.get_sub( ind, {0} );
                                ivg::Matrix obs_epo_nz = obs_epochs.get_sub( ind, {0} );                            
			
                                new_apriori = apr_col_nz.estimate_cpwlf( obs_epo_nz, t_cpwlf );
				
                            }
                            
                            // set rate constraints
                            double rate_cnstr = (double)
                                                params.lookup( param_iter->first )[ j ]["cpwlf"][ "rate_cnstr" ];

                            // modify parameter list
                            for( int d=0; d<t_cpwlf.length(); ++d )
                            {
                                _params.push_back( _params.at( cur_idx.at( i ) ) );
                                _params.back().set_epoch( ivg::Date( t_cpwlf( d ) ) );
                                _params.back().set_rate_cnstr_sigma( rate_cnstr );
                                _params.back().set_apriori( new_apriori(d,0) );

                                if( d == 0 )
                                {
                                    _params.back().set_offset_cnstr_sigma( (double)
                                                                           params.lookup( param_iter->first )[ j ]["polynom"][ "cnstr" ][ 0 ] );
                                }

                                // set pointer to previous (cpwlf) parameter
                                if( d != 0 )
                                {
                                    _params.back().set_offset_cnstr_sigma( 0.0 );
                                }
                            }
                            
                            // delete "original" offset parameter from param_list and apriori matrix
			     
			    _params.erase( _params.begin()+cur_idx.at( i ) );
                            apriori.rem_c( cur_idx.at(i) );
                        }
                    }
                }               
            }
        }
    }
    
    if(gmm && setup["PARAMS"].exists("merge_twin_zwd") && (bool)setup["PARAMS"]["merge_twin_zwd"] ){
        merge_zwd_param( trf, *gmm );
    }
    if(gmm && setup["PARAMS"].exists("merge_twin_grad") && (bool)setup["PARAMS"]["merge_twin_grad"] ){
        merge_grad_param( trf, *gmm );
    }
    
#if DEBUG_VLBI >=2
   cerr << "--- void Param_list::modify_parameterization(Setting &, ivg::Trf &, ivg::Crf &, ivg::Ls_solution & )" << " : " << tim.toc() << " s " << endl;
#endif 
}
// ...........................................................................
void Param_list::set_start_end_epoch(ivg::Date &start_epoch,
                                     ivg::Date &end_epoch)
// ...........................................................................
{
    _start_epoch = start_epoch;
    _end_epoch = end_epoch;
}
// ...........................................................................
bool Param_list::does_include( ivg::paramtype type )
// ...........................................................................
{            
    for(vector<Param>::iterator param_iter = _params.begin(); param_iter < _params.end() ;  param_iter++)
    {
        if(param_iter->get_type() == type)
            return true;
    }
    return false;
}
// ...........................................................................
vector<int> Param_list::get_indexes( vector<ivg::paramtype> types,
                                     string name )
// ...........................................................................
{
    vector<int> out;
    for(int i=0; i<types.size(); i++)
        out.push_back(get_index(types.at(i), name));

    return out;
}
// ...........................................................................
int Param_list::get_index(ivg::paramtype type, string name )
// ...........................................................................
{
    return get_index( Param(type, name) );
}
// ...........................................................................
int Param_list::get_index( ivg::Param parameter )
// ...........................................................................
{
    vector<Param>::iterator param_iter = find( _params.begin(), _params.end(),
                                         parameter );
    if( param_iter != _params.end() )
        return ( param_iter-_params.begin() );
    else
        return -1;
}
// ...........................................................................
std::vector<int> Param_list::get_idx(ivg::paramtype type, std::string name )
// ...........................................................................
{
    std::vector<Param>::iterator param_iter = _params.begin()-1;
    std::vector<int> idx;

    while( param_iter != _params.end() )
    {
        param_iter = find( param_iter+1, _params.end(), Param(type, name) );

        if( param_iter != _params.end() )
            idx.push_back( param_iter - _params.begin() );
    }

    return idx;
}
// ...........................................................................
ivg::Param *Param_list::get_param( int idx )
// ...........................................................................
{
    if(idx == -1 || idx > _params.size()-1)
        throw runtime_error( "ivg::Param Param_list::*get_param( int idx ): Index "+to_string(idx)+" of _params not possible." );
    
    return &(_params.at( idx ));
};

// ...........................................................................
map< string, ivg::Matrix > Param_list::get_station_dependent_param( ivg::paramtype type, param_def pdef, bool prediction, bool sinex )
// ...........................................................................
{
#if DEBUG_VLBI >=2
   cerr << "+++ void Param_list::get_station_dependent_param( ivg::paramtype, param_def )" << endl; 
   tictoc tim;
   tim.tic();
#endif
   
    map< string, vector<double> > delay;
    map< string, vector<double> > std;
    map< string, vector<double> > epo;
    std::vector<double> vec;
    map< string, ivg::Matrix > sta_dep_param;
    
    vector<ivg::Param> parameter;
    if( prediction )
        parameter = _stoch_params;
    else
        parameter = _params;
    
    for( auto &param: parameter )
    {
        if( param.get_type() == type )
        {
            // specify if parameter list was filled based on a SINEX file, since
            // in SINEX definitions the ESTIMATE term referres to the total parameter estimates
            // while the parameter estimates are given as the difference between 
            // the SINEX estimate and the SINEX apriori value.
            if( sinex )
            {
                if( pdef == param_def::aprioris )
                    delay[ param.get_name() ].push_back( param.get_apriori() );
                else if( pdef == param_def::estimates )
                {
                    delay[ param.get_name() ].push_back( param.get_estimate()-param.get_apriori() );            
                    std[ param.get_name() ].push_back( param.get_standard_deviation() );
                }
                else if( pdef == param_def::totals )
                {
                    delay[ param.get_name() ].push_back( param.get_estimate() );
                    std[ param.get_name() ].push_back( param.get_standard_deviation() );
                }    
            }
            else
            {
                if( pdef == param_def::aprioris )
                    delay[ param.get_name() ].push_back( param.get_apriori() );
                else if( pdef == param_def::estimates )
                {
                    delay[ param.get_name() ].push_back( param.get_estimate() );
                    std[ param.get_name() ].push_back( param.get_standard_deviation() );
                }
                else if( pdef == param_def::totals )
                {
                    delay[ param.get_name() ].push_back( param.get_apriori()+param.get_estimate() );            
                    std[ param.get_name() ].push_back( param.get_standard_deviation() );
                }    
            }
                
            epo[ param.get_name() ].push_back( param.get_epoch().get_double_mjd() );
        }
    }
    
    for( auto &it: epo )
    {
        vec = epo[ it.first ];        
        vec.insert( vec.end(), delay[ it.first ].begin(), delay[ it.first ].end() );
                
        if( pdef == param_def::estimates || pdef == param_def::totals )
        {
            vec.insert( vec.end(), std[ it.first ].begin(), std[ it.first ].end() );      
            ivg::Matrix T( vec.begin(), vec.end(), epo[ it.first ].size(), 3 );      
            sta_dep_param[it.first] = T;
        }
        else if( pdef == param_def::aprioris )
        {
            ivg::Matrix T( vec.begin(), vec.end(), epo[ it.first ].size(), 2 );
            sta_dep_param[it.first] = T;
        }            
    }
    
    return sta_dep_param;
    
#if DEBUG_VLBI >=2
   cerr << "--- void Param_list::get_station_dependent_param( ivg::paramtype, param_def ): " 
        << tim.toc() << " s " << endl;
#endif     
};

// ...........................................................................
void Param_list::reduce_params( ivg::Ls_solution & solver)
// ...........................................................................
{
#if DEBUG_VLBI >=2
   cerr << "+++ void Param_list::reduce_params( ivg::Ls_solution &)" << endl; 
   tictoc tim;
   tim.tic();
#endif
   
   // build vector with indexes of params to be reduced end remove params
   // from param-list
   
   std::vector<int> out;
   std::vector<Param>::iterator it = _params.end()-1;
   
   while( it != _params.begin()-1 )
   {
     
      if( it->get_reduce_flag() )
      {
	//if (it->get_typename()=="cbr")
	//{
	//  std::cout << "Reducing: ";
	//  it->show();
	//}
        log<DETAIL>("*** Reducing ") % it->get_name() % " " % it->get_typename() % " " % it->get_order();
	
        out.push_back( it-_params.begin() );
         _params.erase( it );
        
        if(_params.size() == 0)
            throw runtime_error( "void Param_list::reduce_params( ivg::Ls_solution & ): No parameter left in _param_list after reducing.");
      }
      //else if (it->get_typename()=="cbr")
      //{
      //  std::cout << "Not reducing: ";
      //  it->show();
      //}
      

      it--;
   }
   
   // reduce params from LS equation system
   if( out.size() > 0 )
      solver.reduce_params( out );
   
   
#if DEBUG_VLBI >=2
   cerr << "--- void Param_list::reduce_params( ivg::Ls_solution &): " << tim.toc() << " s " << endl;
#endif 
}
// ...........................................................................
void Param_list::insert_param( int idx, ivg::Param param )
// ...........................................................................
{    
    _params.insert( _params.begin()+idx, param );
}
// ...........................................................................
void Param_list::insert_param( std::vector<ivg::Param>::iterator it,
                               ivg::Param param )
// ...........................................................................
{
    // update start and end epoch of the parameter list
    if(param.get_epoch().get_double_mjd() < _start_epoch.get_double_mjd())
        _start_epoch = param.get_epoch();
    else if(param.get_epoch().get_double_mjd() > _end_epoch.get_double_mjd())
        _end_epoch = param.get_epoch();
    
    _params.insert( it, param );
}
// ...........................................................................
void Param_list::remove_param( int idx )
// ...........................................................................
{
    _params.erase( _params.begin()+idx );
}
// ...........................................................................
void Param_list::remove_param( std::vector<ivg::Param>::iterator it )
// ...........................................................................
{
    _params.erase( it );
}

// ...........................................................................
void Param_list::set_estimates( ivg::Matrix x, ivg::Matrix std, bool prediction )
// ...........................................................................
{
    if( prediction )
    {
        for( std::vector<ivg::Param>::iterator it = _stoch_params.begin();
                it != _stoch_params.end(); ++it )
        {
            it->set_estimate( x( it-_stoch_params.begin() ) );
            it->set_standard_deviation( std( it-_stoch_params.begin() ) );
        }        
    }
    else
    {
        for( std::vector<ivg::Param>::iterator it = _params.begin();
                it != _params.end(); ++it )
        {
            it->set_estimate( x( it-_params.begin() ) );
            it->set_standard_deviation( std( it-_params.begin() ) );
        }
    }
}
// ...........................................................................
ivg::Matrix Param_list::extract_apriori_vector()
// ...........................................................................
{
    ivg::Matrix aprioris(_params.size(),1,0.0);
    
    int cnt=0;
    for( std::vector<ivg::Param>::iterator it = _params.begin(); it != _params.end(); ++it )
    {
        aprioris(cnt) = it->get_apriori();
        cnt++;
    }
    return aprioris;
}

// ...........................................................................
ivg::Matrix Param_list::extract_estimate_vector()
// ...........................................................................
{
    ivg::Matrix estimates(_params.size(),1,0.0);
    
    int cnt=0;
    for( std::vector<ivg::Param>::iterator it = _params.begin(); it != _params.end(); ++it )
    {
        estimates(cnt) = it->get_estimate();
        cnt++;
    }
    return estimates;
}

// ...........................................................................
ivg::Matrix Param_list::get_param_data(ivg::paramtype type, int order)
// ...........................................................................
{
    ivg::Matrix data;
    for(std::vector<ivg::Param>::iterator param = _params.begin(); param != _params.end(); ++param)
    {
        if(param->is_type({type}, {order}))
        {
            ivg::Matrix row(1,5,0.0);
            row(0,0) = param->get_epoch().get_double_mjd();
            row(0,1) = param->get_estimate() + param->get_apriori();
            row(0,2) = param->get_estimate();
            row(0,3) = param->get_apriori();
            row(0,4) = param->get_standard_deviation();
            data.append_rows(row);
        } 
    }
    
    return data;
}

// ...........................................................................
ivg::Matrix Param_list::get_param_cpwlf_data( ivg::paramtype type, std::string name )
// ...........................................................................
{
    ivg::Matrix data;
    for(std::vector<ivg::Param>::iterator param = _params.begin(); param != _params.end(); ++param)
    {
        if(param->is_type_name({type}, {name}) && param->is_type({type}, {0}) )
        {
            ivg::Matrix row(1,5,0.0);
            row(0,0) = param->get_epoch().get_double_mjd();
            row(0,1) = param->get_estimate() + param->get_apriori();
            row(0,2) = param->get_estimate();
            row(0,3) = param->get_apriori();
            row(0,4) = param->get_standard_deviation();
            data.append_rows(row);
        } 
    }
    
    return data;
}

// ...........................................................................
ivg::Matrix Param_list::get_param_poly_data( ivg::paramtype type, std::string name )
// ...........................................................................
{
    ivg::Matrix data;
    for(std::vector<ivg::Param>::iterator param = _params.begin(); param != _params.end(); ++param)
    {
        if( param->is_type_name({type}, {name}) && !param->is_type({type}, {0}))
        {
            ivg::Matrix row(1,5,0.0);
            row(0,0) = param->get_epoch().get_double_mjd();
            row(0,1) = param->get_estimate() + param->get_apriori();
            row(0,2) = param->get_estimate();
            row(0,3) = param->get_apriori();
            row(0,4) = param->get_standard_deviation();
            data.append_rows(row);
        } 
    }
    
    return data;
}

// ...........................................................................
bool Param_list::exist_cpwlf_param( ivg::paramtype type, std::string name )
// ...........................................................................
{
    int ctr = 0;
    for(std::vector<ivg::Param>::iterator param = _params.begin(); param != _params.end(); ++param)
    {
        if(param->is_type_name({type}, {name}) && param->is_type({type}, {0}) )
            ctr++;
    }
    
    if( ctr >= 2 )
        return true;
    else
        return false;
}

// ...........................................................................
bool Param_list::exist_polynomial_param( ivg::paramtype type, std::string name )
// ...........................................................................
{
    int ctr = 0;
    bool poly = false;
    for(std::vector<ivg::Param>::iterator param = _params.begin(); param != _params.end(); ++param)
    {
        if(param->is_type_name({type}, {name}) && param->is_type({type}, {0}) )
            ctr++;
        
        if(param->is_type_name({type}, {name}) && !param->is_type({type}, {0}) )
            poly = true;           
    }
    
    if( ctr >= 1 && poly == true )
        return true;
    else
        return false;
}

// ...........................................................................
bool Param_list::exist_stoch_param( ivg::paramtype type, std::string name )
// ...........................................................................
{
    for(std::vector<ivg::Param>::iterator param = _stoch_params.begin(); param != _stoch_params.end(); ++param)
    {        
        if(param->is_type_name({type}, {name}) && param->is_type({type}, {-1}) )
            return true;
    }
    
    return false;   
}

// ...........................................................................
void Param_list::show(string out)
// ...........................................................................
{
    cout << "++++++++++++ Param_list.show(" << out << ") +++++++++++++++++" << endl;
    cout << " _start_epoch: " << setw(8) << left << setprecision(3) << fixed << _start_epoch.get_double_mjd() << endl;
    cout << " MJD / APRIORI+ESTIMATE / APRIORI / ESTIMATE / STD" << endl;
    for(int i=0; i<_params.size(); i++)
    {
        cout << setfill('0') << setw(5) << right << i << " ";
        _params.at(i).show();
    }
    cout << " _end_epoch: " << setw(8) << left << setprecision(3) << fixed << _end_epoch.get_double_mjd() << endl;
    cout << "------------ Param_list.show(" << out << ") -----------------" << endl;
}
// ...........................................................................
void Param_list::show_estimates( bool prediction )
// ...........................................................................
{
    log<RESULT>("***  IDX   #  STATION  TYPE O         MJD                  TOTAL                APRIORI               ESTIMATE     STD-DEVIATION   UNIT");
    
    if( prediction )
    {
        for(int i=0; i<_stoch_params.size(); i++)
        {
            ostringstream ss;
            ss << setfill('0') << setw(4) << right << i << " " << _stoch_params.at(i).get_resultline(true);
            log<RESULT>("*** ") % ss.str();
        }        
    }
    else
    {
        for(int i=0; i<_params.size(); i++)
        {
            ostringstream ss;
            ss << setfill('0') << setw(4) << right << i << " " << _params.at(i).get_resultline(true);
            log<RESULT>("*** ") % ss.str();
        }
    }
    log<RESULT>("***  IDX   #  STATION  TYPE O         MJD                  TOTAL                APRIORI               ESTIMATE     STD-DEVIATION   UNIT");
}
// ...........................................................................
void Param_list::create_constraint_conditions( ivg::Ls_solution &solver )
// ...........................................................................
{
#if DEBUG_VLBI >=2
   cerr << "+++ void Param_list::create_constraint_conditions( ivg::Ls_solution & )" << endl; 
   tictoc tim;
   tim.tic();
#endif
   
    // get number of rate constraints: n CPWLF params = n-1 CPWLF rate constraints
    int n = 2*_params.size();
    int nequ = 0;
    
    
    ivg::Matrix B( n, _params.size(), 0.0 );
    ivg::Matrix wgt( n, 1, 0.0 );

    // loop over all params
    for( std::vector<ivg::Param>::iterator it = _params.begin();
            it != _params.end(); ++it )
    {
        // fill B matrix
        double off, rat;
        it->get_cnstr_sigmas( off, rat );
        // offset constraint
        if( off != 0.0 )
        {
            B( nequ,it-_params.begin() ) = 1.0;
            wgt( nequ ) = 1.0/pow( off,2.0 );

            nequ++;
        }
        // rate constraint => CPWLF
        if( rat != 0.0 )
        {
            // assume that CPWLF params are saved continuously within _params
            if( it != _params.end()-1 && *it == *(it+1) )
            {
                int idx = it-_params.begin();

                ivg::Date d1 = it->get_epoch();
                ivg::Date d2 = (it+1)->get_epoch();
                double dt = ( d2.get_double_mjd()-d1.get_double_mjd() )*24.0;

                B( nequ,idx ) = -1.0;
                B( nequ,idx+1 ) = 1.0;
                wgt( nequ ) = 1.0/pow( dt*rat,2.0 );

                nequ++;
            }
        }
    }

    if( nequ > 0 )
    {
        B = B.get_sub( 0,0,nequ-1,_params.size()-1 );
        wgt = wgt.get_sub( 0,0,nequ-1,0 );

        //B.show();
        solver.update_constraints( B,wgt );
        log<INFO>("*** #") % nequ % " constraints added to system due to create_constraint_conditions";
    }
    
    
#if DEBUG_VLBI >=2
   cerr << "--- void Param_list::create_constraint_conditions( ivg::Ls_solution & )" << " : " << tim.toc() << " s " << endl;
#endif  
}
// ...........................................................................
void Param_list::insert_station_velocities( ivg::Ls_neq &neq_solution, ivg::Date &ref_epoch, ivg::Trf &trf )
// ...........................................................................
{
#if DEBUG_VLBI >=2
   cerr << "+++ void Param_list::insert_station_velocities( ivg::Ls_neq &neq_solution, ivg::Date &ref_epoch )" << endl; 
   tictoc tim;
   tim.tic();
#endif 
   
    vector<int> idx_sta_0;
    vector<double> dt_sta;
    for(vector<ivg::Param>::iterator param = begin(); param!= end(); ++param)
    {
        // only if the paramter is a coordinate, order zero
        if(param->is_type({stax,stay,staz},{0}) )
        {
            ivg::Analysis_station * sta_iter;
            trf.get_station(&sta_iter,param->get_name(), ivg::staname::description);
            
            ivg::Matrix vel = sta_iter->get_vel(ref_epoch);
            double apriori = vel((int)param->get_type()) / ivg::param_unit_fac.at(param->get_type());
            
            idx_sta_0.push_back(param - begin());
            dt_sta.push_back((param->get_epoch().get_double_mjd() - ref_epoch.get_double_mjd())/365.25);           
	    ivg::Matrix xyz0=sta_iter->get_xyz(param->get_epoch());
	   
	   
	    param->set_apriori(param->get_apriori()+apriori*(ref_epoch.get_double_mjd()-param->get_epoch().get_double_mjd())/365.25);
            param->set_epoch(ref_epoch);

	    
	   
            
            param = _params.insert( param+1, ivg::Param(param->get_type(), param->get_name(), ref_epoch, apriori, 1));
	    
        }
    }
    
    //adjust normal equation system
    neq_solution.trafo_params2linear(idx_sta_0, dt_sta);
   
#if DEBUG_VLBI >=2
   cerr << "--- void Param_list::insert_station_velocities( ivg::Ls_neq &neq_solution, ivg::Date &ref_epoch )" << " : " << tim.toc() << " s " << endl;
#endif   
}
// ...........................................................................
void Param_list::transform_cpwlf2offsetrate( ivg::Ls_neq &neq_solution, ivg::Date &t_0 )
// ...........................................................................
{
#if DEBUG_VLBI >=2
   cerr << "+++ void Param_list::transform_cpwlf2offsetrate( ivg::Ls_neq &neq_solution, ivg::Date &t_0 )" << endl; 
   tictoc tim;
   tim.tic();
#endif 
    
    vector<ivg::paramtype> eop_types = {xpo,ypo,ut1,nutx,nuty};
    for(auto &eop: eop_types)
    {
        vector<int> idx = get_idx(eop,"EOP");
        // in case of e.g. XPO , two parameter with order = 0 have to exist
        if(idx.size() == 2 && get_param(idx.at(0))->get_order() == 0 && get_param(idx.at(1))->get_order() == 0)
        {
            log<INFO>("*** Running cpwlf to offset/rate transformation for ") % get_param(idx.at(0))->get_name() % " " % get_param(idx.at(0))->get_typename();

            // get both affected parameter
            ivg::Param *p1 = get_param(idx.at(0));
            ivg::Param *p2 = get_param(idx.at(1));
            
            // epoch to be transformed on
            double t0 = t_0.get_double_mjd();
            // interval epochs left and right
            vector<double> t_i = {p1->get_epoch().get_double_mjd(), p2->get_epoch().get_double_mjd()};

            // adjusting parameter_list 
            // setting new apriori for offset parameter (p1)
            double dt_21 = t_i.at(1)-t_i.at(0);
            double offset_apri = ((t_i.at(1)-t0)/dt_21) * p1->get_apriori() + ((t0-t_i.at(0))/dt_21) * p2->get_apriori();
            p1->set_apriori(offset_apri);
            // setting new parameter epoch
            p1->set_epoch(t_0);

            // in case of XPO, YPO, UT we have to transform to offset + rate
            if(p1->is_type({xpo,ypo,ut1},{0}))
            {
                // transformation of the neq system
                // idx: positions of parameter in the neq_system
                // t_i: interval-epochs of the parameter
                // t0: epoch to be transformed on
                neq_solution.trafo_cpwlf2other(ivg::trafoto::offsetrate, idx, t_i, t0);

                // calculate rate between both sampling points
                double rate_apri = (-1.0/dt_21) * p1->get_apriori() + (1.0/dt_21) * p2->get_apriori();

                double fac_offset = offset_apri*ivg::param_unit_fac.at(p1->get_type());
                double fac_rate = rate_apri*ivg::param_unit_fac.at(p1->get_type());
                log<DETAIL>("*** cpwlf2offsetrate with dt: ") %  dt_21 % " | Offset: " % fac_offset % " | Rate: " % fac_rate ;
                // adjustment of parameter for rate
                p2->set_apriori(rate_apri);
                p2->set_order(1);

                // set ne epoch for p2
                p2->set_epoch(t_0);
            }
            // in case of NUTX, NUTY we have to transform to const (offset only)
            else if(p1->is_type({nutx,nuty},{0}))
            {
                neq_solution.trafo_cpwlf2other(ivg::trafoto::offset, idx, t_i, t0);
                remove_param(idx.at(1));
            }
        }
    }   
   
#if DEBUG_VLBI >=2
   cerr << "--- void Param_list::transform_cpwlf2offsetrate( ivg::Ls_neq &neq_solution, ivg::Date &t_0 )" << " : " << tim.toc() << " s " << endl;
#endif   
}
// ...........................................................................
void Param_list::create_nnr_nnt_equations( Setting &setup, ivg::Trf &trf,
        ivg::Crf &crf, ivg::Ls_solution &solver )
// ...........................................................................
{
#if DEBUG_VLBI >=2
   cerr << "+++ void Param_list::create_nnr_nnt_equations( Setting &, ivg::Trf &,ivg::Crf &, ivg::Ls_solution & )" << endl; 
   tictoc tim;
   tim.tic();
#endif
   
    Setting &params = setup[ "no_net_cnstr" ];
    
    int sta_order = 0;
    if((bool)params["velocities"]["apply"] == true)
        sta_order = 1;

    map<string, vector<ivg::paramtype> > assignment =
    {
        { "stations", {ivg::paramtype::stax, ivg::paramtype::stay, ivg::paramtype::staz}},
        { "sources", {ivg::paramtype::ra,ivg::paramtype::dec}}
    };

    // create list of        ...
    vector<string>sta_names = trf.get_station_names( ivg::staname::ivs_name );
    vector<string>src_names = crf.get_source_names( ivg::srcname::iers );
    vector<string>names;

    // loop over all necessary config files entries in no_net-block
    for( map<string, vector<ivg::paramtype> >::iterator param_iter = assignment.begin(); param_iter != assignment.end(); ++param_iter)
    {
        if( (bool)params.lookup( param_iter->first )["apply"] )
        {

            // select relevant param names
            if( param_iter->first == "stations" )
                names = sta_names;
            else if(param_iter->first == "sources")
                names = src_names;

	   
	    
            // loop over list of this config-file entry
            string selected_group_name = params.lookup( param_iter->first )[0].getName();

            // group index 0  -> all stations, source or what ever
            //             >0 -> this group
            //             <0 -> all but this group
            //int name_idx = params[ param_iter->first ][ selected_group_name ];
            int name_idx = params.lookup( param_iter->first ).lookup( selected_group_name );
	    
            vector<string> param_lst;

            // select all param names of this class (stations, sources, ... ) within this session
            if( name_idx == 0)
            {
                param_lst = names;
       
            }
            // select only those param-names within the requested group
            else if( name_idx > 0 )
            {
                //for( int i=0;i<setup["groups"][selected_group_name][ name_idx-1 ].getLength(); ++i )
                for( int i=0;i<setup["groups"].lookup(selected_group_name)[ name_idx-1 ].getLength(); ++i )
                    //param_lst.push_back( setup["groups"][selected_group_name][ name_idx-1 ][ i ] );
                    param_lst.push_back( setup["groups"].lookup(selected_group_name)[ name_idx-1 ][ i ] );
            }
            // remove entries of given group from the entire list
            else if( name_idx < 0)
            {
                name_idx *= -1;
                param_lst = names;
                //for( int i=0; i<setup["groups"][selected_group_name][ name_idx-1].getLength();
                for( int i=0; i<setup["groups"].lookup(selected_group_name)[ name_idx-1].getLength();
                        ++i )
                {
                    vector<string>::iterator iter = find( param_lst.begin(), param_lst.end(),
                                                          //(const char *)setup["groups"][selected_group_name][ name_idx-1][ i ] );
                                                          (const char *)setup["groups"].lookup(selected_group_name)[ name_idx-1][ i ] );
                    if (iter != param_lst.end())
                        param_lst.erase( iter );
                }
            }
	    
            if( param_iter->first == "stations" || param_iter->first == "velocities"  )
            {
                
	      ivg::Matrix B( 6*(1+sta_order),_params.size(),0.0 );
	      //ivg::Matrix w( 6*(1+sta_order),1,1.0/pow( (double)params[ param_iter->first ]["sigma"], 2.0 ) );
	      ivg::Matrix w( 6*(1+sta_order),1,1.0/pow( (double)params.lookup( param_iter->first )["sigma"], 2.0 ) );
                
                for(int order=0; order<=sta_order; order++)
                {   
                    int sta_cnt=0;
                    for(int i=0; i<param_lst.size(); i++)
                    {                 
                        vector<int> idx;
                        for(auto &param: _params)
                            if(param.is_type({ivg::paramtype::stax, ivg::paramtype::stay, ivg::paramtype::staz},{order}) && param.get_name() == param_lst.at(i))
                                idx.push_back(get_index(param));

                        idx.erase( unique( idx.begin(), idx.end() ), idx.end() );

			vector<int> idx_o0;
                        for(auto &param: _params)
                            if(param.is_type({ivg::paramtype::stax, ivg::paramtype::stay, ivg::paramtype::staz},{0}) && param.get_name() == param_lst.at(i))
                                idx_o0.push_back(get_index(param));

                        idx_o0.erase( unique( idx_o0.begin(), idx_o0.end() ), idx_o0.end() );
			
                        if( idx.size() == 3 )
                        {        
                           log<DETAIL>("*** NNR / NNT datum station: ") % param_lst.at(i);

                           ivg::Matrix Bi( 6,3,0.0 );
                           ivg::Matrix posvel( 3,1,0.0 );

                           for( int i=0; i<3; ++i )
                           {
                               if( _params.at( idx.at( i ) ).get_type() == ivg::paramtype::stax )
                                   posvel( 0 ) = _params.at( idx_o0.at( i ) ).get_apriori();
                               else if( _params.at( idx.at( i ) ).get_type() == ivg::paramtype::stay )
                                   posvel( 1 ) = _params.at( idx_o0.at( i ) ).get_apriori();
                               else if( _params.at( idx.at( i ) ).get_type() == ivg::paramtype::staz )
                                   posvel( 2 ) = _params.at( idx_o0.at( i ) ).get_apriori();
                           }
                           posvel = posvel/(posvel.norm())(0);
			  
                            // only in case of NNR/NNT on stations, NOT velocities
                           // if(order == 0)
                           // {
                                Bi( 0,0 ) = 1.0;
                                Bi( 1,1 ) = 1.0;
                                Bi( 2,2 ) = 1.0;
                                Bi( 4,0 ) = -posvel( 2 );
				//}
                            Bi( 5,0 ) =  posvel( 1 );
                            Bi( 3,1 ) =  posvel( 2 );
                            Bi( 5,1 ) = -posvel( 0 );
                            Bi( 3,2 ) = -posvel( 1 );
                            Bi( 4,2 ) =  posvel( 0 );
			    int offset=6*order;
                            B.set_sub( offset,idx.at(0),Bi( ":",0 ) );
                            B.set_sub( offset,idx.at(1),Bi( ":",1 ) );
                            B.set_sub( offset,idx.at(2),Bi( ":",2 ) );

                            sta_cnt++;
                        }
                        else if( idx.size() == 1 && idx.at(0) != -1 )
                            log<WARNING>("!!! Unexpected length of parameters for NNR/NNT condition: ") % idx.size() % " for " % param_lst.at(i) % " instead of 3";
                    }
                   
                    if(order == 0)
                        log<INFO>("*** #6 constraints added to system due to NNR/NNT on ") % sta_cnt % " of " % sta_names.size() % " station positions";
                    else if(order == 1)
                        log<INFO>("*** #6 constraints added to system due to NNR/NNT on ") % sta_cnt % " of " % sta_names.size() % " station velocities";
                }
		
                solver.update_nnt_nnr_constraints( B,w );
                
            }

            if( param_iter->first == "sources" )
            {
                ivg::Matrix B( 3,_params.size(),0.0 );
                //ivg::Matrix w( 3,1,1.0/pow( (double)params[ param_iter->first ]["sigma"], 2.0 ) );
                ivg::Matrix w( 3,1,1.0/pow( (double)params.lookup( param_iter->first )["sigma"], 2.0 ) );

                int src_cnt=0;
                ivg::Matrix pos( 3,1,0.0 );
                for(int i=0; i<param_lst.size(); i++)
                {
                    vector<int> idx = get_indexes(assignment[ param_iter->first ], param_lst.at(i) );
                    //copy ( idx.begin () , idx.end () , ostream_iterator <short >( cout , ", " ) ); cout << endl;
                    idx.erase( unique( idx.begin(), idx.end() ), idx.end() );
                    //copy ( idx.begin () , idx.end () , ostream_iterator <short >( cout , ", " ) ); cout << endl;

                    double ra, dec;
                    if( idx.size() == 2 )
                    {
		       ivg::Source *sour;
		       crf.get_source(&sour, (string) param_lst.at(i));
		       
                       
                       if(sour->get_n_obs_in_sess()>=(int)params["sources"]["min_obs"])
			 {
			   log<DETAIL>("*** NNR datum source: ") % param_lst.at(i) % " " % sour->get_n_obs_in_sess() % " " % sour->get_n_sessions();
			 ivg::Matrix Bi( 3,2,0.0 );

                       for( int i=0; i<2; ++i )
                       {
                          if( _params.at( idx.at( i ) ).get_type() == ivg::paramtype::dec )
                             dec = _params.at( idx.at( i ) ).get_apriori();
                          else if( _params.at( idx.at( i ) ).get_type() == ivg::paramtype::ra )
                             ra  = _params.at( idx.at( i ) ).get_apriori();
                       }
                    
//                       cerr << " RA: " << ra << " / DEC: " << dec << endl;
                       // RA constraints
                       Bi( 0,0 ) = -sin( dec )*cos( dec )*cos( ra );
                       Bi( 1,0 ) = -sin( dec )*cos( dec )*sin( ra );
                       Bi( 2,0 ) = pow( cos( dec ) ,2.0 );
                       
                       // DEC constraints
                       Bi( 0,1 ) =  sin( ra );
                       Bi( 1,1 ) = -cos( ra );
                       Bi( 2,1 ) = 0.0;
                                              
                       B.set_sub( 0,idx.at(0),Bi( ":",0 ) );
                       B.set_sub( 0,idx.at(1),Bi( ":",1 ) );
                       
                       src_cnt++;
			 }
                    }
                    else if(idx.size() == 1 && idx.at(0) == -1)
                    {
                        // it's ok if a defining source is not found within the session
                    }
                    else
                        throw runtime_error("void Param_list::create_nnr_nnt_equations( ... ): Unexpected length of parameters for NNR condition");
                }
                
                solver.update_nnt_nnr_constraints( B,w );
                log<INFO>("*** #3 constraints added to system due to NNR on ") % src_cnt % " of " % src_names.size() % " sources";
            }
        }
    }
    
#if DEBUG_VLBI >=2
   cerr << "--- void Param_list::create_nnr_nnt_equations( Setting &, ivg::Trf &,ivg::Crf &, ivg::Ls_solution & )" << " : " << tim.toc() << " s " << endl;
#endif 
}
// ...........................................................................
ivg::Matrix Param_list::generate_station_nnr_nnt_B(vector<string> &stations, ivg::Matrix &w)
// ...........................................................................
{
    ivg::Matrix B( 6,_params.size(),0.0 );

    for(int i=0; i<stations.size(); i++)
    {
        vector<int> idx = get_indexes({ivg::paramtype::stax, ivg::paramtype::stay, ivg::paramtype::staz}, stations.at(i) );
        //copy ( idx.begin () , idx.end () , ostream_iterator <short >( cout , ", " ) ); cout << endl;
        idx.erase( unique( idx.begin(), idx.end() ), idx.end() );
        //copy ( idx.begin () , idx.end () , ostream_iterator <short >( cout , ", " ) ); cout << endl;

        if( idx.size() == 3 )
        {                        
           ivg::Matrix Bi( 6,3,0.0 );
           ivg::Matrix pos( 3,1,0.0 );

           for( int i=0; i<3; ++i )
           {
               if( _params.at( idx.at( i ) ).get_type() == ivg::paramtype::stax )
                   pos( 0 ) = _params.at( idx.at( i ) ).get_apriori();
               else if( _params.at( idx.at( i ) ).get_type() == ivg::paramtype::stay )
                   pos( 1 ) = _params.at( idx.at( i ) ).get_apriori();
               else if( _params.at( idx.at( i ) ).get_type() == ivg::paramtype::staz )
                   pos( 2 ) = _params.at( idx.at( i ) ).get_apriori();
           }
           pos = pos/(pos.norm())(0);

           Bi( 0,0 ) = 1.0;
           Bi( 1,1 ) = 1.0;
           Bi( 2,2 ) = 1.0;
           Bi( 4,0 ) = -pos( 2 );
           Bi( 5,0 ) =  pos( 1 );
           Bi( 3,1 ) =  pos( 2 );
           Bi( 5,1 ) = -pos( 0 );
           Bi( 3,2 ) = -pos( 1 );
           Bi( 4,2 ) =  pos( 0 );

           B.set_sub( 0,idx.at(0),Bi( ":",0 ) );
           B.set_sub( 0,idx.at(1),Bi( ":",1 ) );
           B.set_sub( 0,idx.at(2),Bi( ":",2 ) );
        }
    }
    
    return B;
}

// ...........................................................................
void Param_list::merge_zwd_param( ivg::Trf &trf, ivg::Lsa &solver ){
// ...........................................................................
    
#if DEBUG_VLBI >=2
   cerr << "+++ void Param_list::merge_zwd_param( ivg::Trf &trf, ivg::Lsa &solver )" << endl; 
   tictoc tim;
   tim.tic();
#endif
    
    ivg::Matrix* A0 = solver.get_design0_ptr();
    ivg::Matrix* A  = solver.get_design_ptr();
   
    std::vector<int> remove;
    
    std::vector< std::vector<ivg::Analysis_station*> > twin_list = trf.get_station_twin_list();
    // loop over each observatory
    for(std::vector<ivg::Analysis_station*> twins: twin_list ){
        
        if( twins.size() > 1 ){
        
//            std::string newName = twins[0]->get_name(ivg::staname::lettercode);

            // All zwd params are merged into the columns of the first station in the twin list
            // The columns of other station at the observatory are delted afterwards
            std::vector<int> zwd_sta = get_idx( ivg::paramtype::zwd, twins[0]->get_name(ivg::staname::ivs_name) );
            for(unsigned i = 1; i < twins.size(); ++i){

//                newName += twins[i]->get_name(ivg::staname::lettercode);
                std::vector<int> zwd_sta_other = get_idx(ivg::paramtype::zwd, twins[i]->get_name(ivg::staname::ivs_name) );

                log<INFO>("*** merging zwd params of ") % twins[0]->get_name(ivg::staname::ivs_name) % " and " % twins[i]->get_name(ivg::staname::ivs_name);
                
                // merge each zwd (cpwlf) parameter at both stations 
                for(int t = 0; t < zwd_sta.size(); ++t){
                    remove.push_back(zwd_sta_other[t]);

                    // merge columns of jacobian containing zwd param of both stations to one column 
                    for(int r = 0; r < A0->rows(); ++r){
                        if( (*A0)(r, zwd_sta[t]) == 0.0){
                            (*A0)(r, zwd_sta[t]) = (*A0)(r, zwd_sta_other[t]) ;
                        }
                        if( (*A)(r, zwd_sta[t]) == 0.0){
                            (*A)(r, zwd_sta[t]) = (*A)(r, zwd_sta_other[t]) ;
                        } 
                    }
                }
            } 

//            // rename params
//            for(int p :  zwd_sta){
//                _params[p].set_name(newName);
//            }
        
        }
        
    }
        
    // delete params from param_list and jacobian
    std::sort(remove.begin(), remove.end() );
    int c(0);
    for(int& i: remove){
        //std::cout << "remove " <<_params[i-c].get_name() << std::endl;
        this->remove_param(i-c);
        ++c;
    }
    
    std::vector<int> empty;
    solver.remove_data(empty, remove);
    
#if DEBUG_VLBI >=2
   cerr << "--- void Param_list::merge_zwd_param( ivg::Trf &trf, ivg::Lsa &solver )" << " : " << tim.toc() << " s " << endl;
#endif 
    
}


// ...........................................................................
void Param_list::merge_grad_param( ivg::Trf &trf, ivg::Lsa &solver ){
// ...........................................................................
    
#if DEBUG_VLBI >=2
   cerr << "+++ void Param_list::merge_zwd_param( ivg::Trf &trf, ivg::Lsa &solver )" << endl; 
   tictoc tim;
   tim.tic();
#endif
    
    ivg::Matrix* A0 = solver.get_design0_ptr();
    ivg::Matrix* A  = solver.get_design_ptr();
   
    std::vector<int> remove;
    
    std::vector< std::vector<ivg::Analysis_station*> > twin_list = trf.get_station_twin_list();
    // loop over each observatory
    for(std::vector<ivg::Analysis_station*> twins: twin_list ){
        
        if( twins.size() > 1 ){
        
//            std::string newName = twins[0]->get_name(ivg::staname::lettercode);

            // All grad params are merged into the columns of the first station in the twin list
            // The columns of other station at the observatory are delted afterwards
            std::vector<int> egr_sta = get_idx( ivg::paramtype::egr, twins[0]->get_name(ivg::staname::ivs_name) );
	    std::vector<int> ngr_sta =get_idx( ivg::paramtype::ngr, twins[0]->get_name(ivg::staname::ivs_name) );
            for(unsigned i = 1; i < twins.size(); ++i){

//                newName += twins[i]->get_name(ivg::staname::lettercode);
                std::vector<int> egr_sta_other = get_idx(ivg::paramtype::egr, twins[i]->get_name(ivg::staname::ivs_name) );
		std::vector<int> ngr_sta_other = get_idx(ivg::paramtype::ngr, twins[i]->get_name(ivg::staname::ivs_name) );
		
                log<INFO>("*** merging grad params of ") % twins[0]->get_name(ivg::staname::ivs_name) % " and " % twins[i]->get_name(ivg::staname::ivs_name);
                
                // merge each zwd (cpwlf) parameter at both stations 
                for(int t = 0; t < egr_sta.size(); ++t){
                    remove.push_back(egr_sta_other[t]);
		    remove.push_back(ngr_sta_other[t]);
		    
                    // merge columns of jacobian containing zwd param of both stations to one column 
                    for(int r = 0; r < A0->rows(); ++r){
                        if( (*A0)(r, egr_sta[t]) == 0.0){
                            (*A0)(r, egr_sta[t]) = (*A0)(r, egr_sta_other[t]) ;
                        }
                        if( (*A)(r, egr_sta[t]) == 0.0){
                            (*A)(r, egr_sta[t]) = (*A)(r, egr_sta_other[t]) ;
                        }
			if( (*A0)(r, ngr_sta[t]) == 0.0){
                            (*A0)(r, ngr_sta[t]) = (*A0)(r, ngr_sta_other[t]) ;
                        }
                        if( (*A)(r, ngr_sta[t]) == 0.0){
                            (*A)(r, ngr_sta[t]) = (*A)(r, ngr_sta_other[t]) ;
                        } 
                    }
                }
            } 

//            // rename params
//            for(int p :  zwd_sta){
//                _params[p].set_name(newName);
//            }
        
        }
        
    }
        
    // delete params from param_list and jacobian
    std::sort(remove.begin(), remove.end() );
    int c(0);
    for(int& i: remove){
        //std::cout << "remove " <<_params[i-c].get_name() << std::endl;
        this->remove_param(i-c);
        ++c;
    }
    
    std::vector<int> empty;
    solver.remove_data(empty, remove);
    
#if DEBUG_VLBI >=2
   cerr << "--- void Param_list::merge_grad_param( ivg::Trf &trf, ivg::Lsa &solver )" << " : " << tim.toc() << " s " << endl;
#endif 
    
}

// ...........................................................................
void Param_list::create_common_clock_equations( Setting &setup, ivg::Ls_solution &solver )
// ...........................................................................
{
#if DEBUG_VLBI >=2
   cerr << "+++ void Param_list::create_nnr_nnt_equations( Setting &, ivg::Trf &,ivg::Crf &, ivg::Ls_solution & )" << endl; 
   tictoc tim;
   tim.tic();
#endif
   
    Setting &params = setup[ "equal_clocks" ];

    if( (bool)params["apply"] )
    {
        std::vector<int> orders;
        for( int i=0; i < params["orders"].getLength(); ++i )
            orders.push_back( params["orders"][i] );

        ivg::Matrix B( orders.size()*params["stations"].getLength(),_params.size(),0.0 );
        ivg::Matrix w( orders.size()*params["stations"].getLength(),1, 1.0/pow( (double)params["sigma"],2.0) );
        int cnt = 0;
        std::string sta1, sta2; 
        
        for( int j=0; j < params["stations"].getLength(); ++j )
        {      
            sta1 = (const char *)params["stations"][j][0];
            sta2 = (const char *)params["stations"][j][1];
                                    
            if( !( exist_polynomial_param(ivg::paramtype::clo,sta1) && exist_polynomial_param(ivg::paramtype::clo,sta2)) )
                    throw runtime_error( "!!! An equal clock cannot be estimated for stations "+sta1+" and "+sta2+
                                         " is not possible. Maybe, one of these stations is selected as reference clock." );                

            std::vector<int> idx_sta1, idx_sta2;
            for( std::vector<ivg::Param>::iterator param = begin(); param != end(); ++param )
            {
                if( param->is_type( {ivg::paramtype::clo},orders ) && param->get_name() == sta1 )
                    idx_sta1.push_back( param - begin() );      
            }
            for( std::vector<ivg::Param>::iterator param = begin(); param != end(); ++param )
            {
                if( param->is_type( {ivg::paramtype::clo},orders ) && param->get_name() == sta2 )
                    idx_sta2.push_back( param - begin() );
            }
                        
            for( int k=0; k < orders.size(); ++k )
            {
                B(cnt+k,idx_sta1.at(k)) =  1;
                B(cnt+k,idx_sta2.at(k)) = -1;
            }
            
            cnt += orders.size()-1;
        }
        
        solver.update_constraints( B,w );
        log<INFO>("*** ") % B.rows() % " new constraints added to system due common clock estimation for stations " 
                          % sta1 % " and " % sta2 % " .\n";        
    }

#if DEBUG_VLBI >=2
   cerr << "--- void Param_list::create_nnr_nnt_equations( Setting &, ivg::Trf &,ivg::Crf &, ivg::Ls_solution & )" << " : " << tim.toc() << " s " << endl;
#endif 
}

// ...........................................................................
// ...........................................................................
vector<int> Param_list::_get_group_idx( Setting & groups, Setting & params,
                                        vector<string> names, string selected_group_name, int name_idx,
                                        vector<ivg::paramtype> type )
// ...........................................................................
{

    // group index 0  -> all stations, source or what ever
    //             >0 -> this group
    //             <0 -> all but this group
    vector<string> param_lst;

    // select all param names of this class (stations, sources, ... ) within this session
    if( name_idx == 0)
    {
        param_lst = names;
    }

    // select only those param-names within the requested group
    else if( name_idx > 0 )
    {
        for( int i=0; i<groups.getLength(); ++i )
            param_lst.push_back( groups[ i ] );
    }

    // remove entries of given group from the entire list
    else if( name_idx < 0)
    {
        name_idx *= -1;
        param_lst = names;
        for( int i=0; i<groups.getLength(); ++i )
        {
            vector<string>::iterator iter = find( param_lst.begin(), param_lst.end(),
                                                  (const char *)groups[ i ] );
            if (iter != param_lst.end())
                param_lst.erase( iter );
        }
    }

    // create list of indices which are affected by cnf-entry
    vector<int> cur_idx;
    for(int i=0; i<param_lst.size(); i++)
    {
        vector<int> fixis = get_indexes( type,  param_lst.at(i) );
        cur_idx.insert( cur_idx.end(), fixis.begin(), fixis.end() );

        // sort index vector and remove duplicates as well as -1
        sort( cur_idx.begin(), cur_idx.end() );
        cur_idx.erase( unique( cur_idx.begin(), cur_idx.end() ), cur_idx.end() );
        if( cur_idx.at( 0 ) == -1 )
            cur_idx.erase( cur_idx.begin() );
    }

    return cur_idx;
}

// ...........................................................................
void Param_list::_insert_breaks( ivg::Lsa * solver, ivg::paramtype type )
// insert clock breaks based on offset of clock polynomial => no CPWLF or higher
// order polynomials must have been inserted prior to this method
// ...........................................................................
{
    ifstream inStream;
    
    // the map "_breaks[type]" is initialzed either with vgosdb or ngs
    if(!_breaks[type].empty())
    {
        if( type == ivg::paramtype::cbr )
            log<INFO>("*** [#") % _breaks[type].size() % "] clock breaks existent";
        else if( type == ivg::paramtype::atbr )
            log<INFO>("*** [#") % _breaks[type].size() % "] CPWLF breaks for atmospheric parameters existent";
                
        for(auto &cb: _breaks[type])
        {
            std::string station = cb.first;
            double jd = cb.second;

	    if (jd <2400000.5)
	      jd+=2400000.5;
           // try to find indexes of clock parameters for this station, there
           // shoould be excatly one. If none is found (idx.size() == 0),
           // the station might have been eliminated. If more than one is found
           // higher order polynomials or CPWLF might have been introduced
           // already. In the first case, we simply skip this station, in the latter
           // case we abort the solution.
            vector<int> idx;
            if( type == ivg::paramtype::cbr ){
	      idx = get_idx( ivg::paramtype::clo,station );
	    }
            else if( type == ivg::paramtype::atbr )
                idx = get_idx( ivg::paramtype::zwd, station );
            
           if( idx.size() != 0 )  
           {
              // if there is something different than one single clo parameter, throw error
              if( idx.size() > 1 )
                 throw runtime_error( "void Param_list::_insert_breaks( ivg::Lsa * solver, ivg::paramtype type ): Wrong number of clock offsets in solution for station "+station );

              int index = idx.at( 0 );
      
              vector<int> idx_br = get_idx( type, station );
	      if( idx_br.size() > 0 )
                index = idx_br.back();
              // insert break to equation system
              if( type == ivg::paramtype::cbr )
                log<INFO>("*** Insert clock break: ") % station % " " % ( jd-2400000.5 );
              else if( type == ivg::paramtype::atbr )
                log<INFO>("*** Insert CPWLF break for zwd parameter: ") % station % " " % ( jd-2400000.5 );
	      
	      solver->insert_break( index,jd-2400000.5  );
	  
              // add new parameter based on original offset
              ivg::Param br = _params.at( idx.at( 0 ) );

              br.set_type( type );
              br.set_epoch( ivg::Date( jd-2400000.5) );
	    
	      //std::cout << br.get_resultline(false) <<endl;
              _params.push_back( br );
           }
        }
    }
    else
        log<INFO>("*** No clock breaks or CPWLF breaks for atmospheric parameter (ZWDs) existent.");
    
    inStream.close();
}

}

