#include "lsc.h"
#include "param.h"

namespace ivg
{

// ===========================================================================
// 		constructors and destructor
// ===========================================================================

// ...........................................................................
Lsc::Lsc()
// ...........................................................................
{
}

// ...........................................................................
Lsc::Lsc( ivg::Lsa * lsa, ivg::Matrix B )
// ...........................................................................
{
#ifdef DEBUG_SOLVER
    cerr << " +++ Lsc::Lsc( ivg::Lsa * )" << endl;
    tictoc tim;
    tim.tic();
#endif

    _lsa = lsa;
    _B0.eye( _lsa->get_nobs() );
    _Qy.eye( _lsa->get_nobs() );
    _Qw.eye( _lsa->get_nobs() );    
    
#ifdef DEBUG_SOLVER
    cerr << "--- Lsc::Lsc( ivg::Lsa * ) "
         << ": " << tim.toc() << " s " << endl;
#endif
}

// ...........................................................................
Lsc::~Lsc()
// ...........................................................................
{
}


// ===========================================================================
//                              METHODS
// ===========================================================================

// ...........................................................................
void Lsc::find_n_stoch_params()
// ...........................................................................
{
//#ifdef DEBUG_SOLVER
//    cerr << "+++ void Lsc::find_n_stoch_params()" << endl;
//    tictoc tim;
//    tim.tic();
//#endif
    
    for( int i=0; i<_G.size(); ++i )
    {
        ivg::Matrix G = _G.at(i);
//        std::vector<int> idx = G.transpose().sum_col().find_idx( eq, 0.0 );
//        G.rem_r(idx);     
        _nsparam.push_back( G.cols() );
    }
        
//#ifdef DEBUG_SOLVER
//    cerr << "--- void find::get_n_stoch_params()"
//         <<": "<<tim.toc()<<" s "<<endl;
//#endif
}


// ...........................................................................
void Lsc::solve( )
// ...........................................................................
{
#ifdef DEBUG_SOLVER
    cerr << "+++ void Lsc::solve( )" << endl;
    tictoc tim;
    tim.tic();
#endif
                
    ivg::Matrix A = *_lsa->get_design0_ptr();
    ivg::Matrix oc = *_lsa->get_obs0_ptr();  
    
    // modify the weights of the ivg::Lsa object and solve the NEQ to
    // estimate the deterministic parameters
    // (1) calculate the variance-covariance matrix
    ivg::Matrix S = _B0* _Qy *_B0.transpose() + _Qw;    

    // (2) calculate the weight matrix
    ivg::Matrix W;
    W.eye( S.rows() );
    W.solve_neq( S );
    
     // (3) calculate least squares solution vector via ivg::Lsa object
    _lsa->set_weights( W );
    _lsa->reinit();
    _lsa->solve_neq( 0,1.0e-9,3.0,false,0.0 );
    _x = _lsa->get_parameters();
    _Sxx = _lsa->get_vcm() / _lsa->calc_posterior_vfac();
    
      // (3b) alternativ 'classical' least squares approach
//    ivg::Matrix N = A.transpose()*W*A;
//    ivg::Matrix xs = A.transpose()*W*oc;
//    xs.solve_neq(N);
      
    // (4) caluclate residuals of determinic-only residuals (in ivg::Lsa object)
    // [_v0 are scaled by Cholesky-decomposition (i.e., after full decorrelation) 
    //  and needed for data snooping;
    //  _v are the re-scaled residuals needed to the solve the LSC problem]
    _v = _lsa->get_resid();
    _v0 = _lsa->_r;
    
    // (5) estimate the stochastic parameters    
    _y = _Qy* _B0.transpose()* W* ( _v );  
    
    // (5b) calculate the projection matrix P and the covariance matrix of 
    // the estimated stochastic parameters    
    ivg::Matrix P = A* _Sxx* A.transpose()* W;
    ivg::Matrix E;
    E.eye( S.rows() );
    _Syy = _Qy - _Qy* _B0.transpose()* W* (E-P)* _B0* _Qy;   
    
    // (4b) calculate the post-fit residuals of LSC problem and the corresponding
    // variance covariance matrix; and write them into ivg::Lsa object 
    _v = _Qw* W* ( _v ) * -1.0;  
    _v0 = _Qw* W* ( _v0 ) * -1.0;
      
    _Qv = (_lsa->_A*(_Sxx*-1.0)*_lsa->_A.transpose());
    for(int i=0;i<_Qv.rows();i++)
         _Qv(i,i)+=1.0;    
    
    _lsa->set_resid( _v0 );
    _lsa->set_CFM_resid( _Qv );  
    
#ifdef DEBUG_SOLVER
    cerr << "--- void Lsc::solve( )"
         <<": "<<tim.toc()<<" s "<<endl;
#endif
}

// ...........................................................................
void Lsc::calc_stoch_design( std::vector<int> idx )
// ...........................................................................
{
#ifdef DEBUG_SOLVER
    cerr << "+++ void Lsc::calc_stoch_design( )" << endl;
    tictoc tim;
    tim.tic();
#endif  
    
    _B0.resize(0,0,0.0);   
    
    for( int i = 0; i < idx.size(); ++i )
    {
        ivg::Matrix tmp = _lsa->get_design0_ptr()->get_col( idx.at(i) );
        tmp.to_diag_matrix();

//        cerr << "\n *************************** \n B: " << endl;
//        tmp.cout_size();
//        cerr << "rank: " << tmp.rank() << endl;
//
//        cerr << "\n BG: " << endl;
//        (tmp * _G.at(i)).cout_size();
//        cerr << "rank: " << (tmp * _G.at(i)).rank() << endl;
        
        _B0.append_cols( tmp * _G.at(i) );   
        
//        _lsa->get_design_ptr()->rem_c(idx.at(i)-i);
//        _lsa->get_design0_ptr()->rem_c(idx.at(i)-i);
    }        
    
#ifdef DEBUG_SOLVER
    cerr << "--- void Lsc::calc_stoch_design( )"
         << ": " << tim.toc() << " s " << endl;
#endif    
}

// ...........................................................................
void Lsc::calc_stoch_VCM( ivg::Matrix Qy )
// ...........................................................................
{
#ifdef DEBUG_SOLVER
    cerr << "+++ void Lsc::calc_stoch_VCM( )" << endl;
    tictoc tim;
    tim.tic();
#endif

    if( _lsa->_full_vcm )
    {
        ivg::Matrix W = *_lsa->get_wgt_ptr();
        ivg::Matrix I;
        I.eye( W.rows() );
        I.solve_neq( W );
        _Qw = I; 
    }
    else
    {    
        // set covariance matrix for random noise as 1 / wgt
        ivg::Matrix W = *_lsa->get_wgt_ptr();
        ivg::Matrix ones(W.rows(),1,1.0);
        _Qw = ones.div_elem(W);
        _Qw.to_diag_matrix();
    }
    
    // set covariance matrix for the stochastic parameters y
    _Qy = Qy;
    
    // scale covariance matrix due to numerical reasons
//    _scale = pow( sqrt(2)*1e-2 / ivg::c, 2.0 );
    
//    _scale = pow( 1e-2 / ivg::c, 2.0 );
    
//    _scale = pow( sqrt(2*1e-2) / ivg::c, 2.0 );
//    _scale = sqrt(2)*1e-2 / ivg::c;
//    _scale = 1.0 / pow( ivg::param_unit_fac.at( ivg::paramtype::zwd ),2.0 );
//    _scale = pow( 1.0 / ivg::c* 1e-1 ,2.0 );
//    _scale = 1e-21;
//    _Qy *= _scale;
        
    
    _Qy.get_sub(0,0,9,9).show();
    _Qw.get_sub(0,0,9,9).show();
    
    
#ifdef DEBUG_SOLVER
    cerr << "--- void Lsc::calc_stoch_VCM( )"
         <<": "<<tim.toc()<<" s "<<endl;
#endif
}

// ...........................................................................
void Lsc::set_G( std::vector< ivg::Matrix > G )
// ...........................................................................
{
#ifdef DEBUG_SOLVER
    cerr << "+++ void Lsc::set_G( std::vector< ivg::Matrix > G )" << endl;
    tictoc tim;
    tim.tic();
#endif   

    _G = G;
    find_n_stoch_params();
    
#ifdef DEBUG_SOLVER
    cerr << "--- void Lsc::set_G( std::vector< ivg::Matrix > G )" << tim.toc() << " s " << endl;
#endif
}

// ...........................................................................
void Lsc::update_constraints(const ivg::Matrix & B_new,
                             const ivg::Matrix & wgt_new,
                             const ivg::Matrix r_new)
// ...........................................................................
{
#ifdef DEBUG_SOLVER
    cerr << "+++ void Lsc::update_constraints( const ivg::Matrix &"
            << "const ivg::Matrix &,const ivg::Matrix & )" << endl;
    tictoc tim;
    tim.tic();
#endif   

    _lsa->update_constraints( B_new, wgt_new, r_new );
    
#ifdef DEBUG_SOLVER
    cerr << "--- void Lsc::update_constraints( const ivg::Matrix & B_new,const "
         << "ivg::Matrix & wgt_new ): "<<tim.toc()<<" s " << endl;
#endif
}

// ...........................................................................
void Lsc::reduce_params(std::vector<int> idx)
// ...........................................................................
{
#ifdef DEBUG_SOLVER
    cerr << "+++ void Lsc::reduce_param( std::vector<int> ) " << endl;
    tictoc tim;
    tim.tic();
#endif
    
    _lsa->reduce_params( idx );
    
#ifdef DEBUG_SOLVER
    cerr << "--- void Lsc::reduce_param( std::vector<int> )"
         << ": "<<tim.toc()<<" s " << endl;
#endif
}

// ...........................................................................
void Lsc::fix_param(std::vector<int> & idx)
// ...........................................................................
{
#ifdef DEBUG_SOLVER
    cerr << "+++ void Lsc::fix_param( std::vector<int> & idx ) " << endl;
    tictoc tim;
    tim.tic();
#endif

    _lsa->fix_param( idx );

#ifdef DEBUG_SOLVER
    cerr << "--- void Lsc::fix_param( std::vector<int> & idx ) "
         << ": "<<tim.toc() << " s " << endl;
#endif
}

// ...........................................................................
void Lsc::handle_stoch_param( std::vector<int> & idx )
// ...........................................................................
{
#ifdef DEBUG_SOLVER
    cerr << "+++ void Lsc::handle_stoch_param( std::vector<int> & idx ) " << endl;
    tictoc tim;
    tim.tic();
#endif
   
    // sort indices and delete columns beginning with the largest index
    std::sort(idx.begin(),idx.end());
    
    // calculate stochastic design matrix
    calc_stoch_design( idx );
    
    // delete columns of parameters, which should be estimated stochastically,
    // from the original Jacobian matrix (for deterministic parameters)
//    for( int i = idx.size()-1; i>=0; --i )
//    {       
//        _lsa->get_design_ptr()->rem_c(idx.at(i));
//        _lsa->get_design0_ptr()->rem_c(idx.at(i));
//    }

//    fix_param(idx);
      
        
    _nparam -= idx.size();    

#ifdef DEBUG_SOLVER
    cerr << "--- void Lsc::handle_stoch_param( std::vector<int> & idx ) "
         << ": "<<tim.toc() << " s " << endl;
#endif
}

// ...........................................................................
void Lsc::trafo_params2polynomial(int degree,int idx,double t_new, 
                                  ivg::Matrix & aprioris, ivg::Matrix & new_apriori  )
// ...........................................................................
{
#if DEBUG_SOLVER >=2
    cerr << "+++ void Lsc::trafo_params2polynomial( int degree, int idx, double t_new ) " << endl;
    tictoc tim;
    tim.tic();
#endif
    
    _lsa->trafo_params2polynomial( degree, idx, t_new, aprioris, new_apriori );
    
#if DEBUG_SOLVER >=2
    cerr << "--- void Lsc::trafo_params2polynomial( int degree, int idx, double t_new ) "
         <<": "<< tim.toc()<<" s " << endl;
#endif
}

// ...........................................................................
void Lsc::trafo_params2cpwlf(int idx,ivg::Matrix t_new)
// ...........................................................................
{
#if DEBUG_SOLVER >=2
    cerr <<"+++ void Lsc::trafo_params2cpwlf( int idx, ivg::Matrix t_new ) " << endl;
    tictoc tim;
    tim.tic();
#endif
    
    _lsa->trafo_params2cpwlf( idx, t_new );
    
#if DEBUG_SOLVER >=2
    cerr << "--- void Lsc::trafo_params2cpwlf( int idx, ivg::Matrix t_new )"
         << ": "<<tim.toc()<<" s " << endl;
#endif
}

// ...........................................................................
double Lsc::calc_wrms()
// ...........................................................................
{
#ifdef DEBUG_SOLVER
    cerr << "+++ double Lsc::calc_wrms()" << endl;
    tictoc tim;
    tim.tic();
#endif

    double wrms = _lsa->calc_wrms();

#ifdef DEBUG_SOLVER
    cerr << "--- double Lsc::calc_wrms()"
         << ": " << tim.toc() << " s " << endl;
#endif
    return wrms;
}

// ...........................................................................
double Lsc::calc_rms()
// ...........................................................................
{
#ifdef DEBUG_SOLVER
    cerr << "+++ double Lsc::calc_rms()" << endl;
    tictoc tim;
    tim.tic();
#endif

    double rms = _lsa->calc_rms();
 
#ifdef DEBUG_SOLVER
    cerr << "--- double Lsc::calc_rms()"
         << " : " << tim.toc() << " s " << endl;
#endif
    return rms;
}

// ...........................................................................
double Lsc::calc_posterior_vfac()
// ...........................................................................
{
#ifdef DEBUG_SOLVER
    cerr << "+++ double Lsc::calc_posterior_vfac()" << endl;
    tictoc tim;
    tim.tic();
#endif
    
    double vfac = _lsa->calc_posterior_vfac();
    
#ifdef DEBUG_SOLVER
    cerr << "--- double Lsc::calc_posterior_vfac()"
         << ": " << tim.toc() << " s " << endl;
#endif
    return vfac;
}

// ...........................................................................
void Lsc::data_snooping(double alpha,std::string type, int & noutl, double & perc_outl, double & expect_perc_outliers )
// ...........................................................................
{
#ifdef DEBUG_SOLVER
    cerr<<"+++ void Lsc::data_snooping( double )"<<endl;
    tictoc tim;
    tim.tic();
#endif
    
    // maximal percentage of outliers acceptable (default 25.0%)
    #define MAX_OUT_PERC 50.0
    
    vector<int> cols(_lsa->_A.cols());
    std::iota(std::begin(cols),std::end(cols),0);

    // get degrees of freedom
    double dof = _lsa->_get_degrees_of_freedom();

    // outlier detection for each observation
    double T,q,qtmp,ri,s0,r_max;
    int imax;
    noutl = _lsa->_A_outliers.rows();
    bool test = true;

    perc_outl = 0.0;
    while (test)
    {
        // remove previously detected outliers
        ivg::Matrix r = _v0; 
        if (_lsa->_idx_outliers.size()>0)
           r.set_sub(_lsa->_idx_outliers,{0},ivg::Matrix(_lsa->_idx_outliers.size(),1,0.0));
        
        // find observation with the largest (weighted) residual
        r_max = (r.absD()).max(imax);

        // statistical test
        // case 1: given variance factor - BAARDA test
        if (type=="baarda")
        {
            // calculate test value
            ri = _Qv(imax,imax);
            T = r_max/(1.0*sqrt(ri));
                        
            // calculate quantile (= inverse of cummulative distribution function, 
            // cdf ) of standard normal distribution 
            normal s;
            q = quantile(complement(s,alpha/2));
        }
        // case 2: unknown variance factor - POPE test
        else if (type=="pope")
        {
            // calculate test value
            ri = _Qv(imax,imax);
            s0 = calc_posterior_vfac();
            T = r_max/(s0*sqrt(ri));

            // calculate quantile (= inverse of cummulative distribution function, 
            // cdf ) of tau distribution (using fisher distribution)
            fisher_f f(1,dof);
            qtmp = quantile(complement(f,alpha));
            q = sqrt((dof*qtmp)/(dof-1.0+1.0*qtmp));
        }
        else
            throw runtime_error("Lsa::_data_snooping( double ): select BAARDA test (given variance factor) or POPE test (unknown variance factor)! Exiting!");

        if (T>q)
        {
            noutl++;
            log<INFO>("*** ")%type%": "%T%", q = "%q%
                    " >>> null hypothesis rejected; #outliers: "%noutl;
            
            // calculate percentage of outlier
            perc_outl = (((double)noutl)/((double)_v.rows()))*100.0;
            if(perc_outl > MAX_OUT_PERC+expect_perc_outliers)
                throw runtime_error("Lsa::_data_snooping( double ): More than "+std::to_string(MAX_OUT_PERC+expect_perc_outliers)+"% outliers. No further detection sensible.");
                                  
            _lsa->_idx_outliers.push_back( imax );
            if( _lsa->_idx_outliers.size() == 1 )
            {
                _lsa->_A_outliers = _lsa->_A.get_sub({imax},cols);
                _lsa->_b_outliers = _lsa->_b.get_sub({imax},{0});
                _lsa->_A.set_sub({imax},cols,ivg::Matrix(1,_lsa->_A.cols(),0.0));
                _lsa->_b.set_sub({imax},{0},ivg::Matrix(1,1,0.0));
                _lsa->_A0.set_sub({imax},cols,ivg::Matrix(1,_lsa->_A.cols(),0.0));
                _lsa->_b0.set_sub({imax},{0},ivg::Matrix(1,1,0.0));                  
            }
            else
            {
                _lsa->_A_outliers.append_rows(_lsa->_A.get_sub({imax},cols));
                _lsa->_b_outliers.append_rows(_lsa->_b.get_sub({imax},{0}));
                _lsa->_A.set_sub({imax},cols,ivg::Matrix(1,_lsa->_A.cols(),0.0));
                _lsa->_b.set_sub({imax},{0},ivg::Matrix(1,1,0.0));
                _lsa->_A0.set_sub({imax},cols,ivg::Matrix(1,_lsa->_A.cols(),0.0));
                _lsa->_b0.set_sub({imax},{0},ivg::Matrix(1,1,0.0));
            }       
            // new LSC adjustment
            solve();
        }
        else
        {
            log<INFO>("*** ")%type%": "%T%", q = "%q%
                    " >>> null hypothesis accepted";
            test = false;
        }
    }

    log<INFO>("*** #outliers (data snooping)  ")%noutl;
#ifdef DEBUG_SOLVER
    cerr<<"--- void Lsc::data_snooping( double )"
            <<": "<<tim.toc()<<" s "<<endl;
#endif
}

} // namespace ivg
