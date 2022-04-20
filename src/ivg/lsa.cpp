#include "lsa.h"

namespace ivg
{

// ===========================================================================
// 		constructors and destructor
// ===========================================================================

// ...........................................................................
Lsa::Lsa()
// ...........................................................................
{
}

// ...........................................................................
Lsa::Lsa(const ivg::Matrix& A,const ivg::Matrix b,const ivg::Matrix W,
         const ivg::Matrix tb,const ivg::Matrix tx)
// ...........................................................................
{
#ifdef DEBUG_SOLVER
    cerr<<"+++ Lsa::Lsa( const ivg::Matrix&, const ivg::Matrix, "
            <<"const ivg::Matrix, const ivg::Matrix, const ivg::Matrix ) "<<endl;
    tictoc tim;
    tim.tic();
#endif
    // check dimensions
    int n = A.rows();
    int m = A.cols();

    // b should be a column vector with n entries.
    // W might be a nx1, 1xn vector ar a nxn matrix
    int n_w = W.rows();
    int m_w = W.cols();
    bool w_ok = false;
    if ((n_w==1&&m_w==n)||(n_w==n&&m_w==1))
    {
        if (n_w==1)
            _W = W.transpose();
        else
            _W = W;

        _full_vcm = false;
        w_ok = true;
    }
    if (n_w==n&&m_w==n)
    {
        _full_vcm = true;
        w_ok = true;
        _W = W;
    }

    // x-axis values for observations b: ih none are given we set them 
    // equidistant with a spacing of 1
    if (tb.rows()==0)
    {
        ivg::Matrix tmp(0,1,n-1,1);
        _tb = tmp;
    }
    else
    {
        _tb = tb;
        if (_tb.rows()==1)
            _tb = _tb.transpose();
    }

    // if no x-values for parameters are given, these are set to the average x
    // of the observations
    _tx = tx;
    if (_tx.rows()==0)
        _tx.resize(m,1,_tb.meanD());
    else if (_tx.rows()==1)
        _tx = _tx.transpose();

    // throw runtime error if dimensions do not agree
    if (b.rows()!=n|| !w_ok||_tb.rows()!=n||_tx.rows()!=m)
    {
        stringstream errormessage;
        errormessage<<"Lsa::Lsa( const ivg::Matrix& A, const ivg::Matrix b, "
                <<"const ivg::Matrix W, const ivg::Matrix tb, const "
                <<"ivg::Matrix tx ): dimensions do not agree! "
                <<"A("<<n<<","<<m<<"), "
                <<"b("<<b.rows()<<","<<b.cols()<<"), "
                <<"W("<<n_w<<","<<m_w<<"). "<<" "
                <<"tb("<<tb.rows()<<","<<tb.cols()<<"). "<<" "
                <<"tx("<<tx.rows()<<","<<tx.cols()<<"). "
                <<" Exiting";
        throw runtime_error(errormessage.str());
    }

    // full de-correlation
    if (_full_vcm)
        // cholesky decomposition of W for performing full de-correlation 
        _R.chol(_W);
    else
    {
        // do the same for uncorrelated case
        ivg::Matrix r = _W.pow(0.5); // = 1/sigma_i
        _R = r.diag();
    }

    // and incorporate weights into jacobian and observations => new
    // weight matrix is unity
    _A0 = A;
    _b0 = b;
    _A = _R*_A0;
    _b = _R*_b0;
    
    _btPb = (_b.transpose()*_b)(0);

    _nparam = _A.cols();
    _aplo_vec.resize(_b.rows(),1,0.0);
#ifdef DEBUG_SOLVER
    cerr<<"--- Lsa::Lsa( const ivg::Matrix&, const ivg::Matrix, const ivg::Matrix, const ivg::Matrix, const ivg::Matrix ) "
            <<": "<<tim.toc()<<" s "<<endl;
#endif
}

// ...........................................................................
Lsa::~Lsa()
// ...........................................................................
{
}

// ...........................................................................
void Lsa::resize(int nobs,int nparam,bool full_vcm)
// ...........................................................................
{
#ifdef DEBUG_SOLVER
    cerr<<"+++ void Lsa::resize(int,int,bool):"<<endl;
    tictoc tim;
    tim.tic();
#endif
    _full_vcm = full_vcm;
    
    _A.resize(nobs,nparam,0.0);
    _b.resize(nobs,1,0.0);
    _aplo_vec.resize(nobs,1,0.0);
    _tb.resize(nobs,1,0.0);
    _tx.resize(nparam,1,0.0);
    
    if(_full_vcm)
       _W.resize(nobs,nobs,0.0);
    else
       _W.resize(nobs,1,0.0);
    
    _nparam = nparam;
    _nobs = nobs;
    
#ifdef DEBUG_SOLVER
    cerr<<"--- void Lsa::resize(int,int,bool):"
            <<": "<<tim.toc()<<" s "<<endl;
#endif
}
// ...........................................................................
void Lsa::reinit()
// ...........................................................................
{
#ifdef DEBUG_SOLVER
    cerr<<"+++ void Lsa::reinit():"<<endl;
    tictoc tim;
    tim.tic();
#endif
       // check dimensions
    int n = _A.rows();
    int m = _A.cols();

    
    // W might be a nx1, 1xn vector ar a nxn matrix
    int n_w = _W.rows();
    int m_w = _W.cols();
   
    bool w_ok = false;
    if ((n_w==1&&m_w==n)||(n_w==n&&m_w==1))
    {
        if (n_w==1)
            _W = _W.transpose();

        _full_vcm = false;
    }
    if (n_w==n&&m_w==n)
        _full_vcm = true;

    // x-axis values for observations b: if none are given we set them 
    // equidistant with a spacing of 1
    if (_tb.max()==0.0&&_tb.min()==0.0)
    {
        ivg::Matrix tmp(0,1,n-1,1);
        _tb = tmp;
    }
 
    // if no x-values for parameters are given, these are set to the average x
    // of the observations
    if (_tx.max()==0.0&&_tx.min()==0.0)
        _tx.resize(m,1,_tb.meanD());

    // full de-correlation
    if (_full_vcm)
        // cholesky decomposition of W for performing full de-correlation 
        _R.chol(_W);
    else
    {
        // do the same for uncorrelated case
        ivg::Matrix r = _W.pow(0.5); // = 1/sigma_i
        _R = r.diag();
    }
    
    // and incorporate weights into jacobian and observations => new
    // weight matrix is unity
    // note: Since original Jacobian matrix and observation vector are used for
    //       the Least Squares collocation method, only overwrite 
    //       _A0 and _b0, if they are empty!    
    if( _A0.rows() <= 0 && _A0.cols() <= 0)
    {
        _A0 = _A;
        _b0 = _b;
    }
    
    _A = _R*_A0;
    _b = _R*_b0;
    _aplo_vec.resize(_b.rows(),1,0.0);
    _btPb = (_b.transpose()*_b)(0);

    _nparam = _A.cols();
        
#ifdef DEBUG_SOLVER
    cerr<<"--- void Lsa::reinit():"
            <<": "<<tim.toc()<<" s "<<endl;
#endif
}
// ...........................................................................
void Lsa::build_neq()
// ...........................................................................
{
#ifdef DEBUG_SOLVER
    cerr<<"+++ void Lsa::build_neq():"<<endl;
    tictoc tim;
    tim.tic();
#endif
    // normal equations (including VCM of observations and down-weight 
    // factors for observations)    
    if(_full_vcm)
        _A = _R*_A0; // We need this in case of full VCM
        
    _neq.build_neq(_A,_b,_tx,_nparam);
    
    // apply constraints (weights are already transformed to B)
    if (_B.rows()!=0)
    {
        ivg::Matrix ones(_wB.rows(),1,1.0);
        _neq.set_constraints(_B,ones,_rB);
    }

#ifdef DEBUG_SOLVER
    cerr<<"--- void Lsa::build_neq():"
            <<": "<<tim.toc()<<" s "<<endl;
#endif
}

// ...........................................................................
void Lsa::solve_neq(int outlier_iterations,double threshold,
                    double quantile,bool pre_cond,double lambda)
// ...........................................................................
{
#ifdef DEBUG_SOLVER
    cerr<<"+++ void Lsa::solve_neq( int , double , double )"<<endl;
    tictoc tim;
    tim.tic();
#endif
    tictoc tim;
   tim.tic();
    log<RESULT>("*** Dimension of jacobian matrix: ") % _A.rows() % " x " % _A.cols();
    
    
    // solve iteratively for detecting outliers
    // i = 0: find gross errors only
    // i > 0: check n*sigma as well
    for (int i = 0; i<outlier_iterations+1; ++i)
    {
        // (re)construct normal equations due to servaral reasons
        // * do it for the first time
        // * new outliers
        // * the jacobian matrix might have changed (e.g., parameters have been
        //   fixed or reduced)
      
        build_neq();

        // SOLVE
	
        _neq.solve(ivg::solutiontype::neq_chol,pre_cond,lambda);
	
        _x = _neq.get_parameters();        
        _Sxx = _neq.get_vcm();
        
        // try to detected outliers if requested and not in the final iteration
        if (outlier_iterations>0&&i<outlier_iterations) {
	    _detect_outliers(i+1,threshold,quantile);  
	}
	   
    }

    // calculate residuals (including outliers)
    // restore outliers to check them again
    vector<int> cols(_A.cols());
    std::iota(std::begin(cols),std::end(cols),0);
    if (_idx_outliers.size()>0)
    {
        _A.set_sub(_idx_outliers,cols,_A_outliers);
        _b.set_sub(_idx_outliers,{0},_b_outliers);
    }
    
    _r = (_b-_A*_x); // these are still scaled
    ivg::Matrix Nneq, nneq;
    _neq.get_neq(Nneq,nneq,true);
    _btPb=(_r.transpose()*_r+_x.transpose()*nneq)(0);
    // remove outliers again
    if (_idx_outliers.size()>0)
    {
       _A.set_sub(_idx_outliers,cols,ivg::Matrix(_idx_outliers.size(),_A.cols(),0.0));
       _b.set_sub(_idx_outliers,{0},ivg::Matrix(_idx_outliers.size(),1,0.0));
    }
  
   // calculate cofactor matrix of residuals (_Qvv = I-(_A*_Sxx*_A.transpose()))
   _Qvv = (_A*(_Sxx*-1.0)*_A.transpose());
  
   for(int i=0;i<_Qvv.rows();i++)
   	_Qvv(i,i)+=1.0;
   
    if (_dQvv.rows()!=0)
        _Qvv += _dQvv;
    
    // apply  a posteriori variance factor
    
    _Sxx *= calc_posterior_vfac();
#ifdef DEBUG_SOLVER
    cerr<<"--- void Lsa::solve_neq( int , double , double )"
            <<": "<<tim.toc()<<" s "<<endl;
#endif
}

/* 2015-10-14 -TA- removed; has to be check for campatibility with rest of code
// ...........................................................................
void Lsa::solve(int outlier_iterations,double threshold,double quantile,
                bool pre_cond)
// ...........................................................................
{
#ifdef DEBUG_SOLVER
    cerr<<"+++ void Lsa::solve( int , double , double )"<<endl;
    tictoc tim;
    tim.tic();
#endif

    // solve iteratively for detecting outliers
    // i = 0: find gross errors only
    // i > 0: check n*sigma as well
    for (int i = 0; i<outlier_iterations+1; ++i)
    {
        ivg::Matrix A = _A;
        _x = _b;
        if (_B.rows()!=0)
        {
            A.append_rows(_B);
            //ivg::Matrix e( _B.rows(),1,0.0 );
            _x.append_rows(_rB);
        }
        _x.solve_qr(A,_Sxx);

        // try to detect outliers if requested and not in the final iteration
        if (outlier_iterations>0&&i<outlier_iterations)
            _detect_outliers(i+1,threshold,quantile);
    }

   // calculate residuals (including outliers)
    // restore outliers to check them again
    vector<int> cols(_A.cols());
    std::iota(std::begin(cols),std::end(cols),0);
    if (_idx_outliers.size()>0)
    {
        _A.set_sub(_idx_outliers,cols,_A_outliers);
        _b.set_sub(_idx_outliers,{0},_b_outliers);
    }
    
    _r = (_b-_A*_x); // these are still scaled
    // remove outliers again
    if (_idx_outliers.size()>0)
    {
       _A.set_sub(_idx_outliers,cols,ivg::Matrix(_idx_outliers.size(),_A.cols(),0.0));
       _b.set_sub(_idx_outliers,{0},ivg::Matrix(_idx_outliers.size(),1,0.0));
    }

    // calculate cofactor matrix of residuals
    ivg::Matrix I;
    I.eye(_A.rows());
    _Qvv = I-(_A*_Sxx*_A.transpose());
    if (_dQvv.rows()!=0)
        _Qvv += _dQvv;

    // apply  a posteriori variance factor
    _Sxx *= calc_posterior_vfac();

#ifdef DEBUG_SOLVER
    cerr<<"--- void Lsa::solve( int , double , double )"
            <<": "<<tim.toc()<<" s "<<endl;
#endif
}
*/
// ...........................................................................
void Lsa::data_snooping(double alpha,std::string type, int & noutl, double & perc_outl, double & expect_perc_outliers )
// ...........................................................................
{
#ifdef DEBUG_SOLVER
    cerr<<"+++ void Lsa::data_snooping( double )"<<endl;
    tictoc tim;
    tim.tic();
#endif
    
    // maximal percentage of outliers acceptable (default 25.0%)
    #define MAX_OUT_PERC 25.0
    
    //ivg::Matrix idx(_b.rows(),1,0.0);
    vector<int> cols(_A.cols());
    std::iota(std::begin(cols),std::end(cols),0);

    // get degrees of freedom
    double dof = _get_degrees_of_freedom();

    // outlier detection for each observation
    double T,q,qtmp,ri,s0,r_max;
    int imax;
    noutl = _A_outliers.rows();
    bool test = true;

    perc_outl = 0.0;
    while (test)
    {
        // remove previously detected outliers
        ivg::Matrix r = _r; 
        if (_idx_outliers.size()>0)
           r.set_sub(_idx_outliers,{0},ivg::Matrix(_idx_outliers.size(),1,0.0));
        
        // find observation with the largest (weighted) residual
        r_max = (r.absD()).max(imax);

        // statistical test
        if (type=="baarda")
        {
            // case 1: given variance factor - BAARDA test

            // calculate test value
            ri = _Qvv(imax,imax);
            T = r_max/(1.0*sqrt(ri));

            // calculate quantile (= inverse of cummulative distribution function, 
            // cdf ) of standard normal distribution 
            normal s;
            q = quantile(complement(s,alpha/2));
        }
        else if (type=="pope")
        {
            // case 2: unknown variance factor - POPE test

            // calculate test value
            ri = _Qvv(imax,imax);
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
            perc_outl = (((double)noutl)/((double)_r.rows()))*100.0;
            if(perc_outl > MAX_OUT_PERC+expect_perc_outliers)
                throw runtime_error("Lsa::_data_snooping( double ): More than "+std::to_string(MAX_OUT_PERC+expect_perc_outliers)+"% outliers. No further detection sensible.");
                                  
            _idx_outliers.push_back( imax );
            if( _idx_outliers.size() == 1 )
            {
                _A_outliers = _A.get_sub({imax},cols);
                _b_outliers = _b.get_sub({imax},{0});
                _A.set_sub({imax},cols,ivg::Matrix(1,_A.cols(),0.0));
                _b.set_sub({imax},{0},ivg::Matrix(1,1,0.0));
                _A0.set_sub({imax},cols,ivg::Matrix(1,_A.cols(),0.0));
                _b0.set_sub({imax},{0},ivg::Matrix(1,1,0.0));                  
            }
            else
            {
                _A_outliers.append_rows(_A.get_sub({imax},cols));
                _b_outliers.append_rows(_b.get_sub({imax},{0}));
                _A.set_sub({imax},cols,ivg::Matrix(1,_A.cols(),0.0));
                _b.set_sub({imax},{0},ivg::Matrix(1,1,0.0));
                _A0.set_sub({imax},cols,ivg::Matrix(1,_A.cols(),0.0));
                _b0.set_sub({imax},{0},ivg::Matrix(1,1,0.0));
            }       
            // new adjustment
            solve_neq(0,1);
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
    cerr<<"--- void Lsa::data_snooping( double )"
            <<": "<<tim.toc()<<" s "<<endl;
#endif
}

// ...........................................................................
int Lsa::restore_outliers( double q )
// ...........................................................................
{
#ifdef DEBUG_SOLVER
    cerr<<"+++ void Lsa::restore_outliers( double )"<<endl;
    tictoc tim;
    tim.tic();
#endif

    vector<int> cols(_A.cols());
    std::iota(std::begin(cols),std::end(cols),0);

    ivg::Matrix r = get_resid(); 
    r = r.get_sub(_idx_outliers,{0}); 

    // get standard deviations
    ivg::Matrix std;       
    std = _Qvv.diag();
    std = ( std.get_sub(_idx_outliers,{0}) )^(0.5);
        
    // find residuals larger than sigma_obs *quantile    
    r = r.absD() - std * q;
    std::vector<int> idx_outliers = r.find_idx( gt, 0.0 );
    int noutl = (_idx_outliers.size()-idx_outliers.size());

    // if all outliers are going to be restored, we skip this step, because
    // the same observations would be again detected as outliers in the next step
    if( idx_outliers.size() == 0 )
        noutl = 0;

    log<INFO>("*** #outliers restored  ")%noutl;
        
    if( noutl > 0 )
    {        
        _A.set_sub(_idx_outliers,cols,_A_outliers);
        _b.set_sub(_idx_outliers,{0},_b_outliers);  
        
        for( int j = 0; j<idx_outliers.size(); ++j )
            idx_outliers.at(j) = _idx_outliers.at( idx_outliers.at(j) );
        
        _idx_outliers = idx_outliers;
        
        ivg::Matrix Rinv = _get_chol_wgtmat_inv();
        _A0 = Rinv * _A;
        _b0 = Rinv * _b;  
        _A_outliers = _A.get_sub(_idx_outliers,cols);
        _b_outliers = _b.get_sub(_idx_outliers,{0});
        _A.set_sub(_idx_outliers,cols,ivg::Matrix(_idx_outliers.size(),_A.cols(),0.0));
        _b.set_sub(_idx_outliers,{0},ivg::Matrix(_idx_outliers.size(),1,0.0));
        _A0.set_sub(_idx_outliers,cols,ivg::Matrix(_idx_outliers.size(),_A.cols(),0.0));
        _b0.set_sub(_idx_outliers,{0},ivg::Matrix(_idx_outliers.size(),1,0.0));    

        // new adjustment
        solve_neq(0,1);              
    }   
#ifdef DEBUG_SOLVER
    cerr<<"--- void Lsa::restore_outliers( double )"
            <<": "<<tim.toc()<<" s "<<endl;
#endif
    
    return noutl;
}

// ...........................................................................
void Lsa::_detect_outliers(int iteration,double threshold,
                           double quantile)
// ...........................................................................
{
#ifdef DEBUG_SOLVER
    cerr<<"+++ void Lsa::_detect_outliers( int , double , double )"<<endl;
    tictoc tim;
    tim.tic();
#endif
    
    // restore outliers to check them again
    vector<int> cols(_A.cols());
    std::iota(std::begin(cols),std::end(cols),0);

    if (_idx_outliers.size()>0)
    {
        _A.set_sub(_idx_outliers,cols,_A_outliers);
        _b.set_sub(_idx_outliers,{0},_b_outliers);
        _A_outliers.resize(0,0);
        _b_outliers.resize(0,0);
    }
    _r = (_b-_A*_x); // these are still scaled but contain previously detected outliers

    // first test: all residuals larger than 1 ns
    
    ivg::Matrix r0 = _b0-_A0*_x; //get_resid();
     
    _idx_outliers = r0.abs().find_idx(ge,threshold);
    
    // second test: residuals larger than a multiple of sigma
    if (iteration>1)
    {
        // calculate cofactor matrix of residuals
        ivg::Matrix I;
        I.eye(_A.rows());

        _Qvv = I-(_A*_Sxx*_A.transpose());
        if (_dQvv.rows()!=0)
            _Qvv += _dQvv;

        // standard deviations of residuals
        ivg::Matrix std = (_Qvv.diag()).sqrt();

        // calculate test values => positive if residual is larger than n*sigma
        ivg::Matrix test = _r.abs()-std*quantile;

        // find outliers and add them to outlier list from test #1
        vector<int> out_idx2 = test.find_idx(ge,0.0);
        _idx_outliers.insert(_idx_outliers.end(),out_idx2.begin(),out_idx2.end());
        sort(_idx_outliers.begin(),_idx_outliers.end());

        // remove indixes which appear twice
        _idx_outliers.erase(unique(_idx_outliers.begin(),_idx_outliers.end()),
                            _idx_outliers.end());
    }

    log<INFO>("*** #outlier iteration ")%iteration%": "%_idx_outliers.size();
    
    // zero out jacobian matrix
    if (_idx_outliers.size()>0){
    _A_outliers = _A.get_sub(_idx_outliers,cols);
    _A.set_sub(_idx_outliers,cols,ivg::Matrix(_idx_outliers.size(),_A.cols(),0.0));
    _b_outliers = _b.get_sub(_idx_outliers,{0});
    _b.set_sub(_idx_outliers,{0},ivg::Matrix(_idx_outliers.size(),1,0.0));
    }
#ifdef DEBUG_SOLVER
    cerr<<"--- void Lsa::_detect_outliers( int , double , double )"
            <<": "<<tim.toc()<<" s "<<endl;
#endif
}

// ...........................................................................
double Lsa::_get_degrees_of_freedom()
// ...........................................................................
{
#ifdef DEBUG_SOLVER
    cerr<<"+++ double Lsa::_get_degrees_of_freedom()"<<endl;
    tictoc tim;
    tim.tic();
#endif
    
    // boosted calculation of degree of freedom
     double redundancy;
     for(int i=0; i<_Qvv.rows(); i++)
         if(find(_idx_outliers.begin(), _idx_outliers.end(), i) == _idx_outliers.end())
           redundancy += _Qvv(i,i);

    log<DETAIL>("*** dof: ")%redundancy;

#ifdef DEBUG_SOLVER
    cerr<<"--- double Lsa::_get_degrees_of_freedom()"
            <<": "<<tim.toc()<<" s "<<endl;
#endif
    return redundancy;
}

// ...........................................................................
double Lsa::calc_posterior_vfac()
// ...........................................................................
{
#ifdef DEBUG_SOLVER
    cerr<<"+++ double Lsa::calc_posterior_vfac()"<<endl;
    tictoc tim;
    tim.tic();
#endif
    
    // boosted calculation of posteriori vfac
     double redundancy,fac;
     for(int i=0; i<_Qvv.rows(); i++)
       {
         if(find(_idx_outliers.begin(), _idx_outliers.end(), i) == _idx_outliers.end())
           redundancy += _Qvv(i,i);
       }
     ivg::Matrix r = _r;
     r.rem_r(_idx_outliers);
     fac = (r.transpose()*r)(0);
   
#ifdef DEBUG_SOLVER
    cerr<<"--- double Lsa::calc_posterior_vfac()"
            <<": "<<tim.toc()<<" s "<<endl;
#endif
    return fac/redundancy;
}

// ...........................................................................
ivg::Matrix Lsa::get_resid()
// ...........................................................................
{
#ifdef DEBUG_SOLVER
    cerr<<"+++ ivg::Matrix Lsa::get_resid() "<<endl;
    tictoc tim;
    tim.tic();
#endif
    // reconstruct residuals with original units by in case of a full VCM of
    // the observations, the inverse of the cholesky decomposed weight matrix
    // is used. For a diagonal matrix the residuals are divided by the square
    // roots of the weights.
    ivg::Matrix Rinv = _get_chol_wgtmat_inv();
    ivg::Matrix r = Rinv * _r;
#ifdef DEBUG_SOLVER
    cerr<<"--- ivg::Matrix Lsa::get_resid() "
            <<": "<<tim.toc()<<" s "<<endl;
#endif
    return r;
}

// ...........................................................................
void Lsa::get_resid(ivg::Matrix& v,ivg::Matrix& Qvv)
// ...........................................................................
{
#ifdef DEBUG_SOLVER
    cerr<<"+++ ivg::Matrix Lsa::get_resid( ivg::Matrix,ivg::Matrix ) "<<endl;
    tictoc tim;
    tim.tic();
#endif
    // reconstruct residuals with original units by in case of a full VCM of
    // the observations, the inverse of the cholesky decomposed weight matrix
    // is used. For a diagonal matrix the residuals are divided by the square
    // roots of the weights.
    ivg::Matrix Rinv = _get_chol_wgtmat_inv();
    v = Rinv * _r;
    Qvv = Rinv*_Qvv*Rinv.transpose();
#ifdef DEBUG_SOLVER
    cerr<<"--- ivg::Matrix Lsa::get_resid( ivg::Matrix,ivg::Matrix ) "
            <<": "<<tim.toc()<<" s "<<endl;
#endif
}

// ...........................................................................
void Lsa::get_wgt_matrix(ivg::Matrix& W,ivg::Matrix& factors)
// ...........................................................................
{
    factors = ivg::Matrix(_A.rows(),1,1.0);
    if (_idx_outliers.size()>0)
        factors.set_sub(_idx_outliers,{0}, ivg::Matrix(_idx_outliers.size(),1,1e-20));

    W = _W;
}

// ...........................................................................
ivg::Matrix Lsa::get_weights()
// ...........................................................................
{
    return _W;
}

// ...........................................................................
void Lsa::set_weights(ivg::Matrix W)
// ...........................................................................
{
#ifdef DEBUG_SOLVER
    cerr<<"+++ void Lsa::set_weights( ivg::Matrix )"<<endl;
    tictoc tim;
    tim.tic();
#endif
    // check dimensions
    int n = _A.rows();
    int m = _A.cols();
    ivg::Matrix bold = _b;

    // redo full de-correlation
    if( false )
    {
       ivg::Matrix Rinv_old;
       if (_R.rows()>0)
           Rinv_old = _get_chol_wgtmat_inv();
       else
           Rinv_old.eye(_A.rows());
       _A = Rinv_old*_A;
       _b = Rinv_old*_b;
    }
    else
    {
       _A = _A0;
       _b = _b0;
    }

    // W might be a nx1, 1xn vector ar a nxn matrix
    int n_w = W.rows();
    int m_w = W.cols();
    bool w_ok = false;
    if ((n_w==1&&m_w==n)||(n_w==n&&m_w==1))
    {
        if (n_w==1)
            _W = W.transpose();
        else
            _W = W;

        _full_vcm = false;
        w_ok = true;
    }
    if (n_w==n&&m_w==n)
    {
        _full_vcm = true;
        w_ok = true;
        _W = W;
    }

    if (_full_vcm)
        // cholesky decomposition of W for performing full de-correlation 
        _R.chol(W);
    else
        _R = (W^0.5).diag();

    // full de-correlation
    _A = _R*_A;
    _b = _R*_b;
#ifdef DEBUG_SOLVER
    cerr<<"--- double Lsa::set_weights( ivg::Matrix )"
            <<": "<<tim.toc()<<" s "<<endl;
#endif
}

// ...........................................................................
double Lsa::calc_wrms()
// calculate weighted root means squared post fit residuals
// ...........................................................................
{
#ifdef DEBUG_SOLVER
    cerr<<"+++ double Lsa::calc_wrms()"<<endl;
    tictoc tim;
    tim.tic();
#endif
    double wrms;

    // extract weights 
    ivg::Matrix w;
    if (_full_vcm)
        w = _W.diag();
    else
        w = _W;

    // remove outliers
    vector<int> idx = _get_idx_wo_outliers();
    w = w.get_sub(idx,{0});

    if (_full_vcm)
    {
        // as residuals are still scaled by the weights of the observations due 
        // to full de-correlation, we get the "clean" residuals first and weight 
        // them once again
        ivg::Matrix r = get_resid();
        r = r.get_sub(idx,{0}); // get only non-outliers
        wrms = sqrt((r.transpose()*w.diag()*r)(0)/(w.sum_col())(0));
    }
    else
    {
        ivg::Matrix r = _r.get_sub(idx,{0});
        // only diagonal weights, which are already included in the residuals
        wrms = sqrt((r.transpose()*r)(0)/(w.sum_col())(0));
    }
#ifdef DEBUG_SOLVER
    cerr<<"--- double Lsa::calc_wrms()"
            <<": "<<tim.toc()<<" s "<<endl;
#endif
    return wrms;
}
/*
// TEST: calc_wrms() using Qvv
// ...........................................................................
double Lsa::calc_wrms()
// ...........................................................................
{
#ifdef DEBUG_SOLVER
   cerr << "+++ double Lsa::calc_wrms()" << endl;
   tictoc tim;
   tim.tic();
#endif
    int ncnstr = _B.rows();
    std::vector<int> idx = _wgt_factors.find_idx( eq, 1 );
    ivg::Matrix r;
    double wrms;

    // extract weights
    ivg::Matrix w = _w;
    // apply down-weight factors for outliers
    w = _w.get_sub( idx,{0} );

    ivg::Matrix Rinv, resid, Qvv;
    get_resid( resid, Qvv );
    w = Qvv.diag()^(-1.0);

    r = resid.get_sub( idx,{0} );
    wrms =sqrt( ( r.transpose() * w.diag() *r )(0)/(w.sum_col())(0) );

#ifdef DEBUG_SOLVER
   cerr << "--- double Lsa::calc_wrms()" 
        << ": " << tim.toc() << " s " << endl;
#endif

    return wrms;
}
 */

// ...........................................................................
double Lsa::calc_rms()
// calculate root mean squared post fit residuals w/o outliers
// ...........................................................................
{
#ifdef DEBUG_SOLVER
    cerr<<"+++ double Lsa::calc_rms()"<<endl;
    tictoc tim;
    tim.tic();
#endif
    double rms;

    vector<int> idx = _get_idx_wo_outliers();
    ivg::Matrix r = get_resid();
    r = r.get_sub(idx,{0});
    rms = sqrt((r.transpose()*r)(0)/r.rows());
#ifdef DEBUG_SOLVER
    cerr<<"--- double Lsa::calc_rms()"
            <<": "<<tim.toc()<<" s "<<endl;
#endif
    return rms;
}

// ...........................................................................
void Lsa::fix_param(std::vector<int> & idx)
// ...........................................................................
{
#ifdef DEBUG_SOLVER
    cerr<<"+++ void Lsa::fix_param( std::vector<int> & idx ) "<<endl;
    tictoc tim;
    tim.tic();
#endif
    // sort indices and delete columns beginning with the largest index
    std::sort(idx.begin(),idx.end());
    for (int i = idx.size()-1; i>=0; --i)
    {
        _A.rem_c(idx.at(i));
        _A0.rem_c(idx.at(i));
	_B.rem_c(idx.at(i));
        
    }

    _nparam -= idx.size();
#ifdef DEBUG_SOLVER
    cerr<<"--- void Lsa::fix_param( std::vector<int> & idx ) "
            <<": "<<tim.toc()<<" s "<<endl;
#endif
}


// ...........................................................................
void Lsa::trafo_params2polynomial(int degree,int idx,double t_new, ivg::Matrix & aprioris )
// ...........................................................................
{
#if DEBUG_SOLVER >=2
    cerr<<"+++ void Lsa::trafo_params2polynomial( int degree, int idx, double t_new ) "<<endl;
    tictoc tim;
    tim.tic();
#endif
    // calculate time difference between observations and new parameter epoch
    ivg::Matrix dt = _tb-t_new;

    ivg::Matrix org_col = _A(":",idx);
    ivg::Matrix new_cols(_A.rows(),degree);
            
    // fill new sub-matrix with elementwise product of partial derivative wrt
    // constant parameter times power( time difference, degree )
    for (int i = 0; i<degree; ++i)
    {
        new_cols.set_sub(0,i,org_col.mult_elem(dt^(i+1)));
    }
     
    _A.append_cols(new_cols);

    _nparam += new_cols.cols();

    _A0 = _A;
    _A = _R*_A0;
    
#if DEBUG_SOLVER >=2
    cerr<<"--- void Lsa::trafo_params2polynomial( int degree, int idx, double t_new ) "
            <<": "<<tim.toc()<<" s "<<endl;
#endif
}


// ...........................................................................
void Lsa::trafo_params2polynomial(int degree,int idx,double t_new, 
                                  ivg::Matrix & aprioris, ivg::Matrix & new_apriori  )
// ...........................................................................
{
#if DEBUG_SOLVER >=2
    cerr<<"+++ void Lsa::trafo_params2polynomial( int degree, int idx, double t_new ) "<<endl;
    tictoc tim;
    tim.tic();
#endif
    _A = _A0;
    // calculate time difference between observations and new parameter epoch
    ivg::Matrix dt = _tb-t_new;

    ivg::Matrix A_org_col = _A(":",idx);
    ivg::Matrix A0_org_col = _A0(":",idx);
    ivg::Matrix A_new_cols(_A.rows(),degree);
    ivg::Matrix A0_new_cols(_A0.rows(),degree);
    
    ivg::Matrix apr_col = aprioris( ":",idx );
    ivg::Matrix B( _A.rows(),1,1.0 );
        
    // fill new sub-matrix with elementwise product of partial derivative wrt
    // constant parameter times power( time difference, degree )
    for (int i = 0; i<degree; ++i)
    {
        A_new_cols.set_sub(0,i,A_org_col.mult_elem(dt^(i+1)));
        A0_new_cols.set_sub(0,i,A0_org_col.mult_elem(dt^(i+1)));
        B.append_cols( dt^(i+1) );
    }
  
//    ivg::Matrix x;
    std::vector<int> ind = apr_col.find_idx( ne, 0.0 );
    
    vector<int> ind2( B.cols() );
    std::iota(std::begin(ind2),std::end(ind2),0);
    
    if( ind.begin() != ind.end() )
    {
	
        new_apriori = apr_col.get_sub( ind, {0} );
        
        ivg::Matrix Bn = B.get_sub( ind, ind2 );    
    
        ivg::Matrix S;
        new_apriori.solve_qr( Bn,S );
	
    }
    
    _A.append_cols(A_new_cols);

    _nparam += A_new_cols.cols();

    _A0 = _A;
    _A = _R*_A0;
    
#if DEBUG_SOLVER >=2
    cerr<<"--- void Lsa::trafo_params2polynomial( int degree, int idx, double t_new ) "
            <<": "<<tim.toc()<<" s "<<endl;
#endif
}
// ...........................................................................
void Lsa::trafo_params2cpwlf(int idx,ivg::Matrix t_new)
// ...........................................................................
{
#if DEBUG_SOLVER >=2
    cerr<<"+++ void Lsa::trafo_params2cpwlf( int idx, ivg::Matrix t_new ) "<<endl;
    tictoc tim;
    tim.tic();
#endif
    int pos = _A.cols();
    _A.resize(_A.rows(),_A.cols()+t_new.length(),0.0);
    _A0.resize(_A0.rows(),_A0.cols()+t_new.length(),0.0);

    // find interval to which each observation belongs and calculate factor
    // for prior and following parameter
    for (int i = 0; i<_tb.rows(); i++)
    {
        if (_tb(i)!=0.0)
        {
            int j = 1;
            // increase j until we have the parameter interval which contains
            // the current observation
            while (t_new(j)<_tb(i))
                j++;

            // time differences of obs and preceeding/following paramete, i.e.,
            // inner partial derivative
            double dt0 = t_new(j)-t_new(j-1);
            double dt1 = _tb(i)-t_new(j-1);

            // apply inner partials
            _A(i,pos+j-1) = (1.0-dt1/dt0)*_A(i,idx);
            _A(i,pos+j) = (dt1/dt0)*_A(i,idx);
            _A0(i,pos+j-1) = (1.0-dt1/dt0)*_A0(i,idx);
            _A0(i,pos+j) = (dt1/dt0)*_A0(i,idx);
        }
    }

    // remove old (constant) parameter
    _A.rem_c(idx);
    _A0.rem_c(idx);

    _nparam += (t_new.rows()-1);
    
#if DEBUG_SOLVER >=2
    cerr<<"--- void Lsa::trafo_params2cpwlf( int idx, ivg::Matrix t_new )"
            <<": "<<tim.toc()<<" s "<<endl;
#endif
}

void Lsa::rm_nnt_nnr_constraints()
{
  _B.rem_r(_rows_nnt_nnr);
  _rB.rem_r(_rows_nnt_nnr);
  _wB.rem_r(_rows_nnt_nnr);
  _rows_nnt_nnr={};
  build_neq();
  
}

// ...........................................................................
void Lsa::update_nnt_nnr_constraints(const ivg::Matrix & B_new,
                             const ivg::Matrix & wgt_new,
                             const ivg::Matrix r_new)
// ...........................................................................
{
  int rows1=_B.rows();
  update_constraints( B_new,wgt_new,r_new);
  int rows2=_B.rows();
  for (int i=rows1;i<rows2;i++)
    _rows_nnt_nnr.push_back(i);
  

}

// ...........................................................................
void Lsa::update_constraints(const ivg::Matrix & B_new,
                             const ivg::Matrix & wgt_new,
                             const ivg::Matrix r_new)
// ...........................................................................
{
#ifdef DEBUG_SOLVER
    cerr<<"+++ void Lsa::update_constraints( const ivg::Matrix &"
            <<"const ivg::Matrix &,const ivg::Matrix & )"<<endl;
    tictoc tim;
    tim.tic();
#endif   
    // set right hand side to zero, if omited
    ivg::Matrix rB = r_new;
    if (rB.rows()==0||(rB.rows()==2&&rB.cols()==2))
    {
        rB.resize(B_new.rows(),1);
        rB.zero();
    }
    // check dimensions
    if (_A.cols()!=B_new.cols())
    {
        stringstream errormessage;
        errormessage<<"Lsa::update_constraints( const ivg::Matrix &"
                <<"const ivg::Matrix &,const ivg::Matrix & ) "
                <<"dimensions do not agree! "
                <<"A("<<_A.rows()<<","<<_A.cols()<<"), "
                <<"B_new("<<B_new.rows()<<","<<B_new.cols()<<"). "<<" Exiting";
        throw runtime_error(errormessage.str());
    }
    if (wgt_new.rows()!=B_new.rows())
    {
        stringstream errormessage;
        errormessage<<"Lsa::update_constraints( const ivg::Matrix &"
                <<"const ivg::Matrix &,const ivg::Matrix & ) "
                <<"dimensions do not agree! "
                <<"wgt_new("<<wgt_new.rows()<<","<<wgt_new.cols()<<"), "
                <<"B_new("<<B_new.rows()<<","<<B_new.cols()<<"). "<<" Exiting";
        throw runtime_error(errormessage.str());
    }
    if (rB.rows()!=B_new.rows())
    {
        stringstream errormessage;
        errormessage<<"Lsa::update_constraints( const ivg::Matrix &"
                <<"const ivg::Matrix &,const ivg::Matrix & ) "
                <<"dimensions do not agree! "
                <<"r_new("<<r_new.rows()<<","<<r_new.cols()<<"), "
                <<"B_new("<<B_new.rows()<<","<<B_new.cols()<<"). "<<" Exiting";
        throw runtime_error(errormessage.str());
    }

    // include weights into update of design matrix (full de-correlation)
    // pseudo observation are uncorrelated => no cholesky decomposition of the
    // corrsesponding weight matrix necessary
    ivg::Matrix B = B_new;
    for (int i = 0; i<B_new.rows(); ++i)
    {
        rB(i) *= sqrt(wgt_new(i));
        for (int j = 0; j<B_new.cols(); ++j)
            B(i,j) *= sqrt(wgt_new(i));
    }

    // no pseudo observations until now
    if (_B.rows()==0)
    {
        _B = B;
        _wB = wgt_new;
        _rB = rB;
    }
        // append new constraints to constraint matrices
    else
    {
        _B.append_rows(B);
        _wB.append_rows(wgt_new);
        _rB.append_rows(rB);
    }
#ifdef DEBUG_SOLVER
    cerr<<"--- void Lsa::update_constraints( const ivg::Matrix & B_new,const "
            <<"ivg::Matrix & wgt_new ): "<<tim.toc()<<" s "<<endl;
#endif
}

// ...........................................................................
void Lsa::write_backend_gmm(std::string dir,std::string name,
                            std::vector<int> param_idx)
// ...........................................................................
{
#ifdef DEBUG_SOLVER
    cerr<<"+++ void Lsa::write_backend_gmm( std::string dir, std::string name,"
            <<"std::vector<int> param_idx )"<<endl;
    tictoc tim;
    tim.tic();
#endif
    std::vector<int> rows(_A.rows(),0);

    for (int i = 0; i<_A.rows(); ++i)
        rows.at(i) = i;

    ivg::Matrix A = _A.get_sub(rows,param_idx);


    // reconstruct jacobian with original units by in case of a full VCM of
    // the observations, the inverse of the cholesky decomposed weight matrix
    // is used. For a diagonal matrix the residuals are divided by the square
    // roots of the weights.
    ivg::Matrix Rinv = _get_chol_wgtmat_inv();

    (Rinv*A).save_bin(dir+name+"gmm"+"XY"+"_A.dat");
    (Rinv*_b).save_bin(dir+name+"gmm"+"XY"+"_oc.dat");

    _W.save_bin(dir+name+"gmm"+"XY"+"_wgt.dat");
    
    if (_B.rows()!=0)
    {
        ivg::Matrix B = _B;

        // revoke full de-correlation
        for (int i = 0; i<B.rows(); ++i)
            for (int j = 0; j<B.cols(); ++j)
                B(i,j) /= sqrt(_wB(i));

        std::vector<int> rows2(_B.rows(),0);
        for (int i = 0; i<B.rows(); ++i)
            rows2.at(i) = i;
        B = B.get_sub(rows2,param_idx);

        B.save_bin(dir+name+"gmm"+"XY"+"_cnstr.dat");
        _wB.save_bin(dir+name+"gmm"+"XY"+"_cnstr_wgt.dat");
    }
//    ivg::Matrix out(_idx_outliers);
//    out.save_bin(dir+name+"gmm"+"XY"+"_outliers.dat");

#ifdef DEBUG_SOLVER
    cerr<<"--- void Lsa::write_backend_gmm( std::string dir, std::string name,std::vector<int> param_idx )"
            <<": "<<tim.toc()<<" s "<<endl;
#endif
}

// ...........................................................................
void Lsa::remove_data(std::vector<int> rows,std::vector<int> cols)
// ...........................................................................
{
#ifdef DEBUG_SOLVER
    cerr<<"+++ void Lsa::remove_data( std::vector<int> rows, std::vector<int> cols ) "<<endl;
    tictoc tim;
    tim.tic();
#endif

    _A.rem_r(rows);
    _b.rem_r(rows);
    _A0.rem_r(rows);
    _b0.rem_r(rows);
    _tb.rem_r(rows);
    _aplo_vec.rem_r(rows);
    if (_full_vcm)
    {
        _W.rem_r(rows);
        _W.rem_c(rows);
        _R.chol(_W);
    }
    else
    { 
       _W.rem_r(rows); 
       ivg::Matrix r = _W.pow(0.5); // = 1/sigma_i 
       _R = r.diag(); 
    }

    _A.rem_c(cols);
    _A0.rem_c(cols);
    if (_B.rows()>0)
        _B.rem_c(cols);
#ifdef DEBUG_SOLVER
    cerr<<"--- void Lsa::remove_data( std::vector<int> rows, std::vector<int> cols )"
            <<": "<<tim.toc()<<" s "<<endl;
#endif
}

// ...........................................................................
void Lsa::remove_observations( std::vector<int> rows )
// ...........................................................................
{
#ifdef DEBUG_SOLVER
    cerr << "+++ void Lsa::remove_observations( std::vector<int> rows ) " << endl;
    tictoc tim;
    tim.tic();
#endif

    _A.rem_r(rows);
    _b.rem_r(rows);
    _A0.rem_r(rows);
    _b0.rem_r(rows);
    _tb.rem_r(rows);
    _aplo_vec.rem_r(rows);

    if (_full_vcm)
    {
        _W.rem_r(rows);
        _W.rem_c(rows);
        _R.chol(_W);
    }
    else
    { 
       _W.rem_r(rows); 
       ivg::Matrix r = _W.pow(0.5); // = 1/sigma_i 
       _R = r.diag(); 
    }

#ifdef DEBUG_SOLVER
    cerr << "--- void Lsa::remove_observations( std::vector<int> rows )"
         << ": "<<tim.toc()<<" s " << endl;
#endif
}

// ...........................................................................
void Lsa::reduce_params(std::vector<int> idx)
// reduce prameters according to 
// Niemeier (2002) 'Ausgleichungsrechnung', pp 286-289
// Funcke (1982) 'Verfahren zur Parameterelimination im GMM', AVN., 89, pp 112-122
// ...........................................................................
{
#ifdef DEBUG_SOLVER
    cerr<<"+++ void Lsa::reduce_param( std::vector<int> ) "<<endl;
    tictoc tim;
    tim.tic();
#endif
    
    ivg::Matrix A1 = ivg::Matrix(_A0.rows()+_B.rows(),_A0.cols(),0.0);
    A1.set_sub(0,0,_A0); // no full de-correlation at this point
    ivg::Matrix b = _b0;
    ivg::Matrix ap = _aplo_vec;
    
    // constrain LSA (constraints will be eliminated after reduction)
    //A1.append_rows(_B);
    if (_B.rows()>0) {
     A1.set_sub(_A0.rows(),0,_B);
     
     b.append_rows(_rB);
     
     ap.append_rows(_rB*0);
    }
    
    // weight matrix of the observations plus pseudo observations (which are
    // already fully de-correlated => W_pseudo_obs = unity)

    ivg::Matrix W;
    W.eye(_A.rows()+_B.rows());
    if(_full_vcm)
        W.set_sub(0,0,_W);
    else
        W.set_sub(0,0,_W.diag());
    
    // divide A in two parts: A1 remaines, A2 is reduced
    
    ivg::Matrix A2 = A1(':',idx);
    A1.rem_c(idx);
    
    // method of reduced observation equations (Mueller, 1942)
    // -------------------------------------------------------
    // intermediate matrices
    ivg::Matrix Ps = (A2.transpose()*W*A2);
    
    Ps.inv_scal();
    
    ivg::Matrix U = A2*Ps*(A2.transpose()*W);
   
    // reduced jacobian matrix and vector of observations
    A1 = A1-U*A1;
    b = b-U*b;
   
    ap=ap-U*ap;
    // remove constraint equations
    _A0 = A1.get_sub(0,0,_A.rows()-1,A1.cols()-1);
   
    _b0 = b.get_sub(0,0,_A.rows()-1,0);
    
    if (_B.rows()>0) {
        _B = A1.get_sub(_A.rows(),0,A1.rows()-1,A1.cols()-1);
        _rB = b.get_sub(_A.rows(),0,b.rows()-1,0);
    
    }
    _aplo_vec=ap.get_sub(0,0,_A0.rows()-1,0);
    // calculate fully de-correlated jacobian matrix and observations
    
   
    _A = _R*_A0;
    _b = _R*_b0;
  
    // update of VCM of residuals (see Funke (1982) eq(3-35))
    U=U.get_sub(0,0,_A.rows()-1,_A.rows()-1);
    
    
  
    if (_dQvv.rows()==0)
        _dQvv = U*(-1.0);
    else
        _dQvv += U*(-1.0);
    
    _btPb = (_b.transpose()*_b)(0);
#ifdef DEBUG_SOLVER
    cerr<<"--- void Lsa::reduce_param( std::vector<int> )"
            <<": "<<tim.toc()<<" s "<<endl;
#endif
}

// ...........................................................................
double Lsa::cond_ls()
// ...........................................................................
{
#ifdef DEBUG_SOLVER
    cerr<<"+++ void Lsa::cond_ls() "<<endl;
    tictoc tim;
    tim.tic();
#endif

    // check wether residuals have been calculated
    if (_r.length()==0)
    {
    }

    // apply constraints to jacobian matrix
    ivg::Matrix A = _A;
    A.append_rows(_B);
    // and calculate condition of jacobian matrix
    double k = A.cond();
    log<DETAIL>("*** cond(A): ")%k;

    // determine angle between observations and adjusted observations
    double theta = asin((_r.norm())(0)/(_b.norm())(0));
    log<DETAIL>("*** angle between b and A*x: ")%(theta*180.0/M_PI);

    // calculate condition of least squares problem (see Demmel (1997), Ch. 3.3)
    double k_ls = 2*k/cos(theta)+tan(theta)*pow(k,2);

#ifdef DEBUG_SOLVER
    cerr<<"--- void Lsa::cond_ls()"
            <<": "<<tim.toc()<<" s "<<endl;
#endif
    return k_ls;
}

// ...........................................................................
void Lsa::insert_break(int idx,double epoch)
// ...........................................................................
{
    
    // copy original column
    ivg::Matrix new_col = _A0(":",idx);

    // determine indexes before and after the break
    vector<int> idx1 = _tb.find_idx(lt,epoch);
    vector<int> idx2 = _tb.find_idx(ge,epoch);
    
    // zero out elements before ...
     if (idx2.size()>0)
       {
	 ivg::Matrix z2(idx2.size(),1,0.0);

	 _A0.set_sub(idx2,{idx},z2);
       }
    // ... and after the break
    if (idx1.size()>0)
      {
	ivg::Matrix z1(idx1.size(),1,0.0);
	new_col.set_sub(idx1,{0},z1);
      }
  
    // add new parameter
    _A0.append_cols(new_col);   
    _A = _R*_A0;
    
}

// ...........................................................................
void Lsa::combine_params(vector<int> idx,bool remove)
// ...........................................................................
{
    ivg::Matrix new_col = _A(":",idx.at(0));

    for (int i = 1; i<idx.size(); ++i)
        new_col += _A(":",idx.at(i));

    if (remove)
        _A.rem_c(idx);
}

// ...........................................................................
ivg::Matrix Lsa::_get_chol_wgtmat_inv()
// ...........................................................................
{
    ivg::Matrix Rinv;
    if (_full_vcm)
    {
        Rinv.eye(_b.length());
        Rinv.solve_geeqs(_R);
    }
    else
    {
        Rinv = _W.pow(-0.5);
        Rinv.to_diag_matrix();
    }
    return Rinv;
}

// ...........................................................................
std::vector<int> Lsa::_get_idx_wo_outliers()
// ...........................................................................
{
    ivg::Matrix factors(_A.rows(),1,1.0);
    if (_idx_outliers.size()>0)
        factors.set_sub(_idx_outliers,{0},ivg::Matrix(_idx_outliers.size(),1,1e-20));
        
    std::vector<int> idx = factors.find_idx(eq,1);

    return idx;
}
// ...........................................................................
void Lsa::get_neq(ivg::Matrix & N,ivg::Matrix & n,ivg::Matrix & aplo,bool cnstr)
// ...........................................................................
{
   _neq.get_neq(N,n,cnstr);
  
   
   if ((_W.cols()==1)||(_W.cols()==1))
     aplo=_A0.transpose()*_W.diag()*_aplo_vec;
   else
     aplo=_A0.transpose()*_W*_aplo_vec;
  
} 

// ...........................................................................
void Lsa::get_neq(ivg::Matrix & N,ivg::Matrix & n,bool cnstr)
// ...........................................................................
{
   _neq.get_neq(N,n,cnstr);
}


// ...........................................................................
std::vector<int> Lsa::find_undefined_param()
// ...........................................................................
{
   
  //ivg::Matrix A = ivg::Matrix(_A.rows()+_B.rows(),_A.cols(),0.0);
    //A.set_sub(0,0,_A);
    //if(_B.rows()!=0)
    //  A.set_sub(_A.rows(),0,_B);
    
    //  if(_B.rows()!=0)
    //   A.append_rows(_B);


    std::vector<int> id_rem=_A.absD().sum_col().find_idx(0.0);
    if (_B.rows()>0) {
    ivg::Matrix sumBrow=_B.absD().transpose().sum_col();
 
    for (int i=id_rem.size()-1;i>=0;i--)
      {
	for (int j=0;j<_B.rows();j++)
	  {
	    if ((_B(j,id_rem.at(i))!=0)&&(sumBrow(j)!=fabs(_B(j,id_rem.at(i)))))
	      {
		id_rem.erase(id_rem.begin()+i);
		break;
	      }
	  }
      }
    }
    return id_rem;
    //return A.absD().sum_col().find_idx(0.0);
} 

// ...........................................................................
void Lsa::repalce_param(const ivg::Matrix& x, const ivg::Matrix& Sxx)
// ...........................................................................
{
    if ( x.rows() != _A.cols() )
        log<WARNING> ("!!! repalce_param: number of parameters and columns in A must not change" );
    
    // replace parameters and vcm
    _x = x;
    _Sxx = Sxx;
    
    // calcualte Residuals and vcm of residuals
    _r = (_b-_A*_x); // these are still scaled
    
    // calculate cofactor matrix of residuals (_Qvv = I-(_A*_Sxx*_A.transpose()))
    _Qvv = (_A*(_Sxx*-1.0)*_A.transpose());
    for(int i=0;i<_Qvv.rows();i++)
   	_Qvv(i,i)+=1.0;

    if (_dQvv.rows()!=0)
        _Qvv += _dQvv;
    
    // apply  a posteriori variance factor
    _Sxx *= calc_posterior_vfac();
}


} // namespace ivg
