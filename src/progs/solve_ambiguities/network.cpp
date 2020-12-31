
/* 
 * File:   network.cpp
 * Author: corbin
 * 
 * Created on 12. September 2017, 12:27
 */

#include "network.h"

#if DEBUG
unsigned int Network::count = 0;
#endif

// ----------------------------------------------------------------------------
Network::Network() {
// ----------------------------------------------------------------------------
#if DEBUG
    ID = count;
    std::cerr << "+++ Network::Network() ID: " << ID << std::endl;
    count++;
#endif
}

// ----------------------------------------------------------------------------
Network::Network(ivg::Session* session) {
// ----------------------------------------------------------------------------
#if DEBUG
    ID = count;
    std::cerr << "+++ Network::Network( ivg::Session* session) ID: " << ID << std::endl;
    count++;
#endif
    
    _refSta = string((const char*) (*session->get_setup())["RefClockStationList"]);
    
    // get groupdelay residuals and indices in nc file
    std::map <delaytype, std::vector<Residual> > resids = session->get_residuals();
    for (auto &r : resids[delaytype::group]) {
        if (r.type == residtype::baseline) {
            Baseline tmp;

            // get station names. they are used as key
            std::string::size_type n = r.name.find("-");
            std::string sta1 = r.name.substr(0, n);
            std::string sta2 = r.name.substr(n + 1);

            // get groupdelay residuals and indices in nc file
            tmp.residuals = r.data.get_col(1).get_rows(r.idx_in);
            tmp.obs_idxs_in_ncfile = r.idx_ncfile;

            // reserve memory for matrices that are filled later (with push_back)
            unsigned int m = r.idx_ncfile.size();
            tmp.delays.data_vec_ptr()->reserve(m);
            tmp.delays_sigma.data_vec_ptr()->reserve(m);

            tmp.integerAmbiguities.resize(m, 1, 0.0);
            tmp.shift.resize(m, 1, 0.0);
            tmp.ambiguity_spacing = session->get_ambiguity_spacing().operator()(tmp.obs_idxs_in_ncfile);

            tmp.fixed = false;
            
            // add to map
            _baselines[std::pair<std::string, std::string>(sta1, sta2)] = tmp;

        }
    }

    // get observations
    unsigned counter = 0;
    Setting *setting = session->get_setup();
    for (ivg::Scan & scan : *(session->get_scan_ptr())) {
        for (unsigned i = 0; i < scan.get_nobs(); ++i) {
            std::string sta1, sta2;
            scan.get_obs_ptr(i)->get_station_names(sta1, sta2);

            std::pair<std::string, std::string> key(sta1, sta2);
	    if ((bool)setting->exists("phase_solution")) {
	      if ((*setting)["phase_solution"]){
		_baselines[key].delays.append_rows(scan.get_obs_ptr(i)->get_phase_delay());
		_baselines[key].delays_sigma.append_rows(scan.get_obs_ptr(i)->get_sigma_phase_delay());
	      } else {
		_baselines[key].delays.append_rows(scan.get_obs_ptr(i)->get_group_delay());
		_baselines[key].delays_sigma.append_rows(scan.get_obs_ptr(i)->get_sigma_group_delay());
	      }
	    } else {
	      _baselines[key].delays.append_rows(scan.get_obs_ptr(i)->get_group_delay());
	      _baselines[key].delays_sigma.append_rows(scan.get_obs_ptr(i)->get_sigma_group_delay());
	    }
            ++counter;
        }
    }

    _n_bl = _baselines.size();
    _n_tr = 0;
    
    // get some other information stored in session so that networks works if session is destructed
    _nobs_orig = session->get_nobs_orig();
    _effFreq = session->get_eff_freq();
    _name = session->get_name();
    

}

// ----------------------------------------------------------------------------
Network::~Network() {
// ----------------------------------------------------------------------------
#if DEBUG
    std::cerr << "--- ~Network() ID: " << ID << std::endl;
#endif
}

// ----------------------------------------------------------------------------
void Network::_createTriangles(bool includesRefStation) {
// ----------------------------------------------------------------------------
#if DEBUG
    std::cerr << "... Network::_createTriangles(bool includesRefStation) ID: " << ID << std::endl;
#endif
    
    log<INFO> ("*** creating triangles");
    
    _triangles.clear();
    
    _triangles.reserve(std::round(_n_bl / 2));

    // build adjacency list
    std::map<std::string, std::vector<std::string> > adjacency;

    std::map< std::pair<std::string, std::string >, Baseline>::const_iterator bl;
    for (bl = _baselines.begin(); bl != _baselines.end(); bl++) {
        adjacency[bl->first.first].push_back(bl->first.second);
    }

    // for each station evaluate the adjacencies
    std::map<std::string, std::vector<std::string> >::const_iterator ad;
    for (ad = adjacency.begin(); ad != adjacency.end(); ad++) {
        // station has at least two neighbours
        if (ad->second.size() >= 2) {

            // find neighbours that are adjacent to each other -> Triangle
            for (unsigned short i = 0; i < ad->second.size(); ++i) {
                for (unsigned short j = 0; j < ad->second.size(); ++j) {
                    if (i != j) {
                        if (std::find(adjacency[ad->second[i]].begin(),
                                adjacency[ad->second[i]].end(),
                                ad->second[j]) != adjacency[ad->second[i]].end()) {

                            // ignore found triangle if reference station is not includeed and "includesRefStation" is true
                            if (includesRefStation)
                                if (ad->first.compare(_refSta) != 0 && ad->second[i].compare(_refSta) && ad->second[j].compare(_refSta))
                                    continue;

                            // create Triangle and add it to vector
                            // Note: due to [] access a  new entrie is created if not existent -> may cause error if key is wrong
                            Triangle tri( std::array<Baseline*, 3>({
                                         &_baselines [ std::pair<std::string, std::string> (ad->first,     ad->second[i]) ],
                                         &_baselines [ std::pair<std::string, std::string> (ad->second[i], ad->second[j]) ],
                                         &_baselines [ std::pair<std::string, std::string> (ad->first,     ad->second[j]) ]
                                         }),
                                         std::array<string, 3>({ad->first, ad->second[i], ad->second[j]}));
                            _triangles.push_back(tri);
                        }
                    }
                }

            }

        }
    }
    
    _n_tr = _triangles.size();


}

// ----------------------------------------------------------------------------
void Network::compute_baselinewise_integer_ambiguities(Algorithem alg, ivg::Session* session) {
// ----------------------------------------------------------------------------
#if DEBUG
    std::cerr << "... compute_baselinewise_integer_ambiguities(Algorithem alg) ID: " << ID << std::endl;
#endif

    log<INFO> ("*** computing integer ambiguities for each baseline");
    switch (alg) {
        case(Algorithem::manual):
        {

            log<WARNING> ("*** manual ambiguity resoltion mode");
            // loop over all baselines and get the integer ambiguities saved in session.
            // in statistics ui this matrix is filled manually
            std::map< std::pair<std::string, std::string >, Baseline>::iterator bl;
            for (bl = _baselines.begin(); bl != _baselines.end(); bl++) {
                bl->second.shift = session->get_num_ambig_ptr()->operator()(bl->second.obs_idxs_in_ncfile);
            }

            break;
        }
        case(Algorithem::ahc):
        {
#ifdef USE_AHC

            // loop over all baselines
            std::map< std::pair<std::string, std::string >, Baseline>::iterator bl;
            for (bl = _baselines.begin(); bl != _baselines.end(); bl++) {


                // http://www.alglib.net/translator/man/manual.cpp.html#sub_clusterizersetpoints
                alglib::clusterizerstate s;
                alglib::ahcreport rep;

                alglib::real_2d_array xy; // variable containing the data to be clustred
                unsigned rows = bl->second.residuals.rows();
                xy.setlength(rows, 1);

                for (int i = 0; i < rows; ++i) {
                    xy(i, 0) = bl->second.residuals(i);
                }

                alglib::clusterizercreate(s);
                alglib::clusterizersetpoints(s, xy, 2); // Euclidean distance  (L2 norm), non-squared

                alglib::clusterizersetahcalgo(s, 2); // 2  unweighted average linkage
                alglib::clusterizerrunahc(s, rep); // performs agglomerative hierarchical clustering

                //number of clusters k
                alglib::ae_int_t k = 0;


                //array[NPoints], I-th element contains cluster index  (from 0 to K-1) for I-th point of the dataset.
                alglib::integer_1d_array cidx;
                // array[K]. This array allows  to  convert  cluster  indexes returned by this function to indexes used by  Rep.Z.  J-th
                // cluster returned by this function corresponds to  CZ[J]-th cluster stored in Rep.Z/PZ/PM.
                // It is guaranteed that CZ[I]<CZ[I+1].
                alglib::integer_1d_array cz;
                alglib::clusterizerseparatedbydist(rep, bl->second.ambiguity_spacing.max()*0.75, k, cidx, cz);

                // this map contains for every cluster the corresponding index
                map<int, vector<int>> clusters;
                for (int j = 0; j < cidx.length(); ++j) {
                    clusters[cidx[j]].push_back(j);
                }

                ivg::Matrix clust_mean; // mean of clusters
                int maxp = 0;
                int c_ref = 0; //reference cluster
                //calculate mean and find the refernce cluster
                int count = 0;
                for (auto c : clusters) {
                    //get the mean of the cluster
                    clust_mean.append_rows(bl->second.residuals.get_vec(c.second).meanD());
                    if (maxp < c.second.size()) {
                        maxp = c.second.size();
                        c_ref = c.first;
                    }// if cluster have the same size choose the one closer to zero
                    else if (maxp == c.second.size()) {
                        if (abs(clust_mean(count)) <= abs(clust_mean(c_ref))) {

                            maxp = c.second.size();
                            c_ref = c.first;
                        }
                    }
                    ++count;

                }

                if (k == 1) {
                    log<DETAIL>("only one cluster found -> no ambiguities in baseline ") % bl->first.first % " - " % bl->first.second;
                } else {
                    log<DETAIL>("found ") % k % " clusters. Cluster " % std::to_string(c_ref+1) % " has the most points and is used as reference for baseline " % bl->first.first % " - " % bl->first.second;
                }


                ivg::Matrix d = ivg::Matrix(k, 1, clust_mean(c_ref)) - clust_mean;
                ivg::Matrix a = (d / (bl->second.ambiguity_spacing.max())).round();

                for (unsigned i = 0; i < rows; ++i) {
                    bl->second.shift(i) = a(cidx[i]);
                }

            }

#else
            std::cerr << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
            std::cerr << ">>> AGGLOMERATIVE HIRACHICAL CLUSTERING (AHC) IS NOT SUPPORTED IN THIS VERSION. use most neighbours instead <<<" << std::endl;
            std::cerr << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
#endif 
            break;
        }
        case(Algorithem::mn):
        {

            // loop over all baselines
            std::map< std::pair<std::string, std::string >, Baseline>::iterator bl;
            for (bl = _baselines.begin(); bl != _baselines.end(); bl++) {

                unsigned int rows = bl->second.residuals.rows();

                // threshold for nearest neigbours -> everything within 
                double thres = bl->second.ambiguity_spacing.max() / 40;

                // find for each residual the residuals that are within the threshold -> neighbours
                ivg::Matrix num_neigbours(rows, 1, 0.0);
                for (unsigned int i = 0; i < rows; ++i) {
                    for (unsigned int j = i; j < rows; ++j) {
                        if (i != j) {
                            double dist = bl->second.residuals(i) - bl->second.residuals(j);
                            if (abs(dist) < thres) {
                                ++num_neigbours(i);
                                ++num_neigbours(j);
                            }
                        }

                    }
                }

                // value of residual with most neighbours
                double reference_resid = bl->second.residuals(num_neigbours.abs().maxIdx());

                ivg::Matrix d = ivg::Matrix(rows, 1, reference_resid) - bl->second.residuals;

                // how often is the distance divisible by the amiguity spacing
                ivg::Matrix f = d.div_elem(bl->second.ambiguity_spacing); //float ambigs

                // Breaks float ambigs into an integral and a fractional part. both parts have the same sign as f
                ivg::Matrix fractpart = f.modf(bl->second.shift);
                
                for (unsigned i = 0; i < rows; ++i) {
                    // check weather fractpart is larger than a threshold (necessary due to noise in observations)
                    // if so enlarge integerpart
                    if (abs(fractpart(i)) > 0.5) {
                        if (d(i) >= 0) {
                            bl->second.shift(i) += 1;
                        } else {
                            bl->second.shift(i) -= 1;
                        }
                    }
                }

            }
               
            break;
        }
    }
    
    _shift_residuals_and_delays();
}

// ----------------------------------------------------------------------------
void Network::apply_closure_condition(){
// ----------------------------------------------------------------------------
#if DEBUG
    std::cerr << "... Network::apply_closure_condition() ID: " << ID << std::endl;
#endif
    
    _createTriangles();

    
    if( _n_tr > 0){
        log<INFO> ("*** applying triangle closure condition");
        
        for( Triangle& tri : _triangles){
            tri.close_loop(_refSta);
        }
        
        _shift_residuals_and_delays();
        
    } else {
        log<INFO> ("*** no triangle in network. #triangles ") % _n_tr % " #baselines " % _n_bl;
    }
    
}

// ----------------------------------------------------------------------------
void Network::fill_all_delay_and_integer_matrix( ) {
// ----------------------------------------------------------------------------
#if DEBUG
    std::cerr << "... fill_all_delay_and_integer_matrix ID: " << ID << std::endl;
#endif
    
    _all_integerAmbiguities.resize(_nobs_orig, 666);
    _all_group_delays.resize(_nobs_orig, 1, NAN);
    _all_group_delays_sigmas.resize(_nobs_orig, 1, NAN);
    _all_group_delays_residuals.resize(_nobs_orig, 1, NAN);
       
    std::map< std::pair<std::string, std::string >, Baseline>::const_iterator bl;
    for (bl = _baselines.begin(); bl != _baselines.end(); bl++) {
        for(unsigned int i = 0; i < bl->second.obs_idxs_in_ncfile.size(); ++i){
                _all_integerAmbiguities[ bl->second.obs_idxs_in_ncfile[i] ] = bl->second.integerAmbiguities(i); 
                      _all_group_delays( bl->second.obs_idxs_in_ncfile[i] ) = bl->second.delays(i);
               _all_group_delays_sigmas( bl->second.obs_idxs_in_ncfile[i] ) = bl->second.delays_sigma(i); 
               _all_group_delays_residuals( bl->second.obs_idxs_in_ncfile[i] ) = bl->second.residuals(i); // needed to update residuals in session
        }
    }   
}

// ----------------------------------------------------------------------------
void Network::get_ionospheric_correction(const Network& X_net, const Network& S_net,
                                   ivg::Matrix& delta_tau_x, ivg::Matrix& delta_tau_x_sigma, std::vector<short>& error_flag ){
// ----------------------------------------------------------------------------
#if DEBUG
    std::cerr << "static void Network::get_ionospheric_correction: " << std::endl;
#endif
    
    if( X_net._name.compare( S_net._name ) !=  0){
        log<WARNING> ("!!! X and S net have to be derived from the same session");
    }
    
    unsigned int n = X_net._nobs_orig;
    
    // Matrix containing inosphere correction in the first column and the rate in the second col
    delta_tau_x.resize(n,2, NAN);
    
    // Matrix containing sigmas of correction in first col, and sigma rate in second row
    delta_tau_x_sigma.resize(n,2, NAN);
    
    // ionospheric correction ErrorFlag 0=OK , -1 = Missing, -2 = bad
    error_flag.resize(n, -1);
    
    
    // f = fs^2 / (fx^2 - fs^2)    dimension of f (MHz GHz ...) does not matter as long as fs and fx have the same dimension it is canceled
    ivg::Matrix f(n,1);
    for(unsigned int i = 0; i < n; ++i){
        f(i) = std::pow( S_net._effFreq(i), 2.0) / ( std::pow( X_net._effFreq(i), 2.0) - std::pow( S_net._effFreq(i), 2.0) );
        
        // not all elements in _all_group_delays may be filled due to the vgos useflags
        // therefor set the error flag to 0=OK if there are observations in X and S band
        if( !std::isnan(S_net._all_group_delays(i)) && !std::isnan(X_net._all_group_delays(i)) ){
            // TODO -2 if bad, criterion for bad?
            error_flag[i] = 0; 
        }
    }

    delta_tau_x.set_col(0, (S_net._all_group_delays - X_net._all_group_delays).mult_elem(f) );
    delta_tau_x_sigma.set_col(0, ( (S_net._all_group_delays_sigmas^2.0) + (X_net._all_group_delays_sigmas^2.0) ).sqrt().mult_elem(f) );
    
    // rate is not calculeated yet. It is not used in ascot so far. If you need it fell free and extend the code...
   
}

// ----------------------------------------------------------------------------
void Network::_shift_residuals_and_delays(){
// ----------------------------------------------------------------------------
#if DEBUG
    std::cerr << "... Network::shift_residuals_and_delays() ID: " << ID << std::endl;
#endif

    log<INFO> ("*** shifting groupdelay residuals and groupdelays");
    std::map< std::pair<std::string, std::string >, Baseline>::iterator bl;
    for (bl = _baselines.begin(); bl != _baselines.end(); bl++) {
        // update integer ambiguities
        bl->second.integerAmbiguities += bl->second.shift;
        
        // update delays and residuals
        ivg::Matrix shift = bl->second.shift.mult_elem(bl->second.ambiguity_spacing);  
        bl->second.delays += shift;
        bl->second.residuals += shift;
        
        // after shifting set shift to zero
        bl->second.shift.zero();
    }
    
}

// ----------------------------------------------------------------------------
void Network::update_resids_in_session(ivg::Session* session){
// ----------------------------------------------------------------------------

    std::map <delaytype, std::vector<Residual> > resids = session->get_residuals();
    for (auto &r : resids[delaytype::group]) {
        
           // indices in nc file
           ivg::Matrix inds = r.idx_ncfile;
           
            for(unsigned int i = 0; i < r.data.rows(); ++i){
                // multi-band delays stehen in Spalte 1
                 r.data(i, 1) = _all_group_delays_residuals( inds(i) );
            }

    }
    
    session->set_residuals(resids);
}

// ----------------------------------------------------------------------------
void Network::print_Network() const{
// ----------------------------------------------------------------------------
#if DEBUG
    std::cerr << "... Network::print_Network() const ID: " << ID << std::endl;
#endif
    
    std::cout << "network information -------------------" << std::endl;
    std::cout << "#Baselines: " << _n_bl << std::endl;
    std::map< std::pair<std::string, std::string >, Baseline>::const_iterator it;
    for (it = _baselines.begin(); it != _baselines.end(); it++) {
        std::cout << "  " << it->first.first << " " << it->first.second << std::endl;
        std::cout << "    " << "observations: " << it->second.obs_idxs_in_ncfile.size() << std::endl;
        std::cout << "    " << "mean residual: " << std::fixed  << std::setprecision(1) 
                            <<  it->second.residuals.meanD()*1E+9 << " ns" << std::endl;

        //        std::cout << "  indx in nc file" << std::endl << "  ";
        //        show_vector(it->second.obs_idxs_in_ncfile); 
        //        std::cout << "  groupdelay"  << "  ";
        //        it->second.delays.show();
        //        std::cout << "  groupdelay residuals" << "  ";
        //        it->second.residuals.show();
        //        std::cout << std::endl;
    }


    std::cout << "#Triangles: " << _n_tr << std::endl;
    for (const Triangle& tri : _triangles) {
        tri.print_Triangle();
    }

}



