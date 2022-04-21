#include "ls_neq.h"
#include "ivg_const.h"

namespace ivg
{

// ===========================================================================
// 		constructors and destructor
// ===========================================================================

// ...........................................................................
Ls_neq::Ls_neq()
// ...........................................................................
{
    _sigma0_pre = 1.0;
    _nparam = 0;
}

// ...........................................................................
Ls_neq::Ls_neq(const ivg::Matrix& N,const ivg::Matrix& n,const int nparam,const ivg::Matrix& epochs)
// ...........................................................................
{
#ifdef DEBUG_SOLVER
    cerr<<"+++ Ls_neq::Ls_neq( const ivg::Matrix&, const ivg::Matrix, "
            <<"const ivg::Matrix, const ivg::Matrix, const ivg::Matrix ) "<<endl;
    tictoc tim;
    tim.tic();
#endif

    set_neq(N,n,epochs);

    //    _N.save_bin("/home/iddink/NNR/N_orig.bin");
    //   cerr << "N after reading" << endl;

#ifdef DEBUG_SOLVER
    cerr<<"--- Ls_neq::Ls_neq( const ivg::Matrix&, const ivg::Matrix, const ivg::Matrix ) "
            <<": "<<tim.toc()<<" s "<<endl;
#endif
}

// ...........................................................................
Ls_neq::Ls_neq(const ivg::Ls_neq &other)
// ...........................................................................
{
    _x = other._x;
    _Sxx = other._Sxx;
    _tx = other._tx;
    _nobs = other._nobs;
    _nparam = other._nparam;
    _sigma0_pre = other._sigma0_pre;
    _sigma0_post = other._sigma0_post;
    _btPb = other._btPb;
    _rtPr = other._rtPr;
    _B = other._B;
    _rB = other._rB;
    _wB = other._wB;


    _N = other._N;
    _dN = other._dN;
    _dn = other._dN;
    _n = other._n;
    _scales = other._scales;

}

// ...........................................................................
Ls_neq::~Ls_neq()
// ...........................................................................
{
}

// ===========================================================================
// 		public methods
// ===========================================================================

// ...........................................................................
void Ls_neq::set_neq(const ivg::Matrix& N,const ivg::Matrix& n,const ivg::Matrix& epochs,const int nparam)
// ...........................................................................
{
#ifdef DEBUG_SOLVER
    cerr<<"+++ Ls_neq::set_neq( const ivg::Matrix&, const ivg::Matrix, "
            <<"const ivg::Matrix, const ivg::Matrix, const ivg::Matrix ) "<<endl;
    tictoc tim;
    tim.tic();
#endif

    // throw runtime error if dimensions do not agree
    if (N.rows()!=N.cols()||n.rows()!=N.rows())
    {
        stringstream errormessage;
        errormessage<<"Ls_neq::set_neq( const ivg::Matrix& N, const ivg::Matrix n, "
                <<"const ivg::Matrix tx ): dimensions do not agree! "
                <<"N("<<N.rows()<<","<<N.cols()<<"), "
                <<"n("<<n.rows()<<","<<n.cols()<<"). "<<" Exiting";
        throw runtime_error(errormessage.str());
    }

    _N = N;
    _n = n;

    _tx = epochs;

    _sigma0_pre = 1.0;

    if (nparam!=0)
        _nparam = nparam;
    else
        _nparam = N.rows();

#ifdef DEBUG_SOLVER
    cerr<<"--- Ls_neq::set_neq( const ivg::Matrix&, const ivg::Matrix, const ivg::Matrix ) "
            <<": "<<tim.toc()<<" s "<<endl;
#endif
}

// ...........................................................................
void Ls_neq::build_neq(const ivg::Matrix& A,const ivg::Matrix& b,
                       const ivg::Matrix& epochs,const int nparam)
// ...........................................................................
{
#ifdef DEBUG_SOLVER
    cerr<<"+++ Ls_neq::build_neq( const ivg::Matrix&, const ivg::Matrix, "
            <<"const ivg::Matrix ) "<<endl;
    tictoc tim;
    tim.tic();
#endif

    // throw runtime error if dimensions do not agree
    if (A.rows()!=b.rows())
    {
        stringstream errormessage;
        errormessage<<"Ls_neq::build_neq( const ivg::Matrix& A, const ivg::Matrix b, "
                <<"const ivg::Matrix epochs ): dimensions do not agree! "
                <<"A("<<A.rows()<<","<<A.cols()<<"), "
                <<"b("<<b.rows()<<","<<b.cols()<<"). "<<" Exiting";
        throw runtime_error(errormessage.str());
    }

    // change calculation of AtA (ICLS)
//    _N.is_AtA_of(A);
//    _n.is_Atv(A,b);
    _N = A.transpose()*A;
    _n = A.transpose()*b;
    
    
    _tx = epochs;

    if (nparam!=0)
        _nparam = nparam;
    else
        _nparam = A.cols();

    _sigma0_pre = 1.0;

#ifdef DEBUG_SOLVER
    cerr<<"--- Ls_neq::build_neq( const ivg::Matrix&, const ivg::Matrix, "
            <<"const ivg::Matrix ) : "<<tim.toc()<<" s "<<endl;
#endif
}

// ...........................................................................
void Ls_neq::build_neq(const ivg::Matrix& A,const ivg::Matrix& b,const ivg::Matrix& W,
                       const ivg::Matrix& epochs,const int nparam)
// ...........................................................................
{
#ifdef DEBUG_SOLVER
    cerr<<"+++ Ls_neq::build_neq( const ivg::Matrix&, const ivg::Matrix, "
            <<"const ivg::Matrix, const ivg::Matrix ) "<<endl;
    tictoc tim;
    tim.tic();
#endif
    ivg::Matrix P;
    if (W.rows()==1||W.cols()==1)
        P = W.diag();
    else
        P = W;

    // throw runtime error if dimensions do not agree
    if (A.rows()!=b.rows()||A.rows()!=P.rows()||P.rows()!=P.cols())
    {
        stringstream errormessage;
        errormessage<<"Ls_neq::set_neq( const ivg::Matrix& N, const ivg::Matrix n, "
                <<"const ivg::Matrix tx ): dimensions do not agree! "
                <<"A("<<A.rows()<<","<<A.cols()<<"), "
                <<"W("<<W.rows()<<","<<W.cols()<<"), "
                <<"b("<<b.rows()<<","<<b.cols()<<"). "<<" Exiting";
        throw runtime_error(errormessage.str());
    }

    _N = A.transpose()*P*A;
    _n = A.transpose()*P*b;

    _tx = epochs;

    if (nparam!=0)
        _nparam = nparam;
    else
        _nparam = A.cols();

    _sigma0_pre = 1.0;

#ifdef DEBUG_SOLVER
    cerr<<"--- Ls_neq::build_neq( const ivg::Matrix&, const ivg::Matrix, "
            <<"const ivg::Matrix, const ivg::Matrix ) : "<<tim.toc()<<" s "<<endl;
#endif
}

// ...........................................................................
void Ls_neq::fix_param(std::vector<int> & idx)
// ...........................................................................
{
#ifdef DEBUG_SOLVER
    cerr<<"+++ void Ls_neq::fix_param( std::vector<int> & idx ) "<<endl;
    tictoc tim;
    tim.tic();
#endif

    //   cerr << "[";copy ( idx.begin () , idx.end () , ostream_iterator <int >( cerr , ", " ) ); cerr<< "];" << endl;

    _check_negative_index(idx);

    _N.rem_r(idx);
    _N.rem_c(idx);
    _n.rem_r(idx);
    

    if (_dN.cols()!=0)
    {
        _dN.rem_r(idx);
        _dN.rem_c(idx);
        _dn.rem_r(idx);
    }

    if (_scales.rows()!=0)
    {
        _scales.rem_r(idx);
    }

    _nparam -= idx.size();


#ifdef DEBUG_SOLVER
    cerr<<"--- void Ls_neq::fix_param( std::vector<int> & idx ) "
            <<": "<<tim.toc()<<" s "<<endl;
#endif

}
// ...........................................................................
void Ls_neq::trafo_cpwlf2other(ivg::trafoto type, vector<int> idx, vector<double> t_i , double t_0)
// ...........................................................................
{
#ifdef DEBUG_SOLVER
    cerr<<"+++ vvoid Ls_neq::trafo_cpwlf2other(ivg::trafoto type, vector<int> idx, vector<double> t_i , double t_0)"<<endl;
    tictoc tim;
    tim.tic();
#endif
    
    if (_dN.rows()!=0)
        throw runtime_error("void Ls_neq::trafo_cpwlf2other(ivg::trafoto type, vector<int> idx, vector<double> t_i , double t_0): _dN already set.");
    
    // transformation only possible if two cpwlf are given, not three!
    if(idx.size() == t_i.size() && idx.size() == 2)
    {
        int n_old = idx.size(); // 2 -> 2xcpwlf 
        int n_new = 2; // 2 -> 1xoffset + 1xrate
        double t_1 = t_i.at(0);
        double t_2 = t_i.at(1);
        
        int i = idx.at(0);
        int j = idx.at(1);
        
        ivg::Matrix Ti(n_new, n_old);
        // OLD and seems to be WRONG?!
        //Ti(0,0) = (t_2 - t_0) / (t_2 - t_1);
        //Ti(1,0) = (t_0 - t_1) / (t_2 - t_1);
        //Ti(0,1) = -1.0 / (t_2 - t_1);
        //Ti(1,1) = 1.0 / (t_2 - t_1);
        
        Ti(0,0) = 1.0;
        Ti(1,0) = 1.0;
        Ti(0,1) = (t_1 - t_0);
        Ti(1,1) = (t_2 - t_0);

        ivg::Matrix B; // B [#old_params x #new_params]
        B.eye(_nparam);
        B.set_sub(i,i,Ti);

        if(type == ivg::trafoto::offset)
        {
            B.rem_c(j);
            _nparam--;
        }
        
        _N = B.transpose() * _N * B;
        _n = B.transpose() * _n;
    }
    else
        throw runtime_error("void Ls_neq::trafo_cpwlf2other(ivg::trafoto type, vector<int> idx, vector<double> t_i , double t_0): Transformation only implemented for two sampling points.");
    
    
#ifdef DEBUG_SOLVER
    cerr<<"--- void Ls_neq::trafo_cpwlf2other(ivg::trafoto type, vector<int> idx, vector<double> t_i , double t_0)" <<": "<<tim.toc()<<" s "<<endl;
#endif
}
// ...........................................................................
void Ls_neq::trafo_params2linear(vector<int> idx_0, vector<double> dt)
// ...........................................................................
{
#ifdef DEBUG_SOLVER
    cerr<<"+++ void Ls_neq::trafo_params2linear(vector<int> idx_0, vector<double> dt)"<<endl;
    tictoc tim;
    tim.tic();
#endif
    
    if (_dN.rows()!=0)
        throw runtime_error("void Ls_neq::trafo_params2linear(vector<int> idx_0, vector<double> dt): _dN already set. PROBLEM?!");
    
    // transformation only possible if idx_0 and dt is of same length
    if(idx_0.size() == dt.size())
    {
        int n_old = _N.rows(); // e.g. 10 stations a x,y,z + 8 EOPs + 20 sources a ra,dec
        int n_new = n_old + idx_0.size(); // e.g. idx_0 should be 30 (3*10)
        
        ivg::Matrix T(n_old,n_new,0.0);
             
        // generate transformation matrix T (identity for all param, which will
        // be added by columns that consist of time differences to ref. epoch)
        int c=0;
        for(int r=0; r<T.rows(); r++)
        {
            if(find( idx_0.begin(), idx_0.end(), c ) != idx_0.end())
            {
                T(r,c) = 1.0;
                T(r,c+1) = dt.at(r);
                c++;
            }
            else
                T(r,c) = 1.0;
            
            c++;
        }
        
        _N = T.transpose() * _N * T;
        _n = T.transpose() * _n;
        
        _nparam = _N.rows();
    }
    else
        throw runtime_error("void Ls_neq::trafo_params2linear(vector<int> idx_0, vector<double> dt): Transformation only implemented for two sampling points.");
    
    
#ifdef DEBUG_SOLVER
    cerr<<"--- void Ls_neq::trafo_params2linear(vector<int> idx_0, vector<double> dt)" <<": "<<tim.toc()<<" s "<<endl;
#endif
}
// ...........................................................................
void Ls_neq::trafo_params2polynomial(int deg,int idx,double t_new)
// ...........................................................................
{
}

// ...........................................................................
void Ls_neq::trafo_params2cpwlf(int idx,ivg::Matrix t_new)
// ...........................................................................
{
}
// ...........................................................................
vector<int> Ls_neq::detect_erase_unparameterized()
// ...........................................................................
{
#ifdef DEBUG_SOLVER
    cerr<<"+++ vector<int> Ls_neq::detect_erase_unparameterized()"<<endl;
    tictoc tim;
    tim.tic();
#endif
    
    vector<int> indexes;
    for(int i=0; i<_n.rows(); i++)
        if(_N(i,i) == 0.0) // maybe also _n(i) == 0.0
            indexes.push_back(i);
    
    log<INFO>("*** Erasing ") % indexes.size() % " columns/rows from NEQ due to zeros in NEQ";
    fix_param(indexes);
    
    return indexes;
  
#ifdef DEBUG_SOLVER
    cerr<<"---vector<int>  Ls_neq::detect_erase_unparameterized()" <<": "<<tim.toc()<<" s "<<endl;
#endif
}

// ...........................................................................
void Ls_neq::new_solve()
// ...........................................................................
{
#ifdef DEBUG_SOLVER
    cerr<<"+++ void Ls_neq::new_solve( ) "<<endl;
    tictoc tim;
    tim.tic();
#endif
    cerr<<"MANUELLER ABBRUCH"<<endl;
    abort();

    if (_N.rows()==0||_N.cols()==0)
    {
        stringstream errormessage;
        errormessage<<"Ls_neq::solve( ivg::solutiontype type, bool pre_cond ): _N.size() = [0,0] ";
        throw runtime_error(errormessage.str());
    }

    //   cerr << "CONDITION: " << _N.cond() << endl;

    // apply constraints; no impact on right hand side as the constraint 
    // equation equal zero
    ivg::Matrix N = _N;
    ivg::Matrix n = _n;
    if (_dN.rows()!=0)
    {
        N += _dN;
        n += _dn;
    }


    ivg::Matrix inverse;
    inverse = N;
    inverse.inv();

    _x = inverse * n;
    _Sxx = inverse;

    if (_scales.rows()!=0)
        _scale_system(_Sxx,_x,_scales);

    //    N.save_bin("/home/iddink/NNR/N_final.bin");
    //    _x.save_bin("/home/iddink/NNR/n_final.bin");


#ifdef DEBUG_SOLVER
    cerr<<"--- void Ls_neq::new_solve( ) "
            <<": "<<tim.toc()<<" s "<<endl;
#endif
}

// ...........................................................................
void Ls_neq::solve(ivg::solutiontype type,bool pre_cond,double lambda)
// ...........................................................................
{
#ifdef DEBUG_SOLVER
    cerr<<"+++ void Ls_neq::solve( type, bool ) "<<endl;
    tictoc tim;
    tim.tic();
#endif
    
    if (_N.rows()==0||_N.cols()==0)
    {
        stringstream errormessage;
        errormessage<<"Ls_neq::solve( ivg::solutiontype type, bool pre_cond ): _N.size() = [0,0] ";
        throw runtime_error(errormessage.str());
    }
    
    // apply constraints; no impact on right hand side as the constraint 
    // equation equal zero   
    ivg::Matrix N = _N;
    ivg::Matrix n = _n;
    
    if (_dN.rows()!=0)
    {
        if (_scales.rows()!=0)
        {
            log<INFO>("*** Scaling _dN");
            ivg::Matrix trash = _scales;
            _scale_system(_dN,_dn,_scales);
        }

        N += _dN;
        n += _dn;
    }
    
    // tikhonov regularization
    if( lambda != 0.0)
    {
        ivg::Matrix E; E.eye(N.rows()); E = E*pow(lambda,2);
        N += E;
    }

    // sort indices and delete columns beginning with the largest index
    // pre-conditioning
    ivg::Matrix S;
    if (pre_cond)
    {
        S = N.diag();
        S = S^(-0.5);
        S = S.diag();
    }
    else
        S.eye(_N.rows());

    // scale normal equation matrix so that the main diagonal consists
    // of ones, if pre_cond = true
    //   double cond1 = N.cond();
    
    N = S*N*S;
    //   double cond2 = N.cond();
    //   log<INFO>("*** condition(N): ") % cond1 % "; condition(S*N*S) " % cond2;
    
    // solve normal equations
    _x = S*n;
    
    if (type==ivg::solutiontype::neq_chol)
        _x.solve_neq(N);
    if (type==ivg::solutiontype::neq_lu)
        _x.solve_neq(N);
    
    // re-scale estimated parameters to original unity (changed by
    // pre-conditioning)
    _x = _x.mult_elem(S.diag());

    // calculate covariance matrix of estimated parameters
    _Sxx.eye(N.rows());
    _Sxx.solve_neq(N);
    _Sxx = S*_Sxx*S;

    // if the system has been scaled with _scales, re-scale now!
    if (_scales.rows()!=0)
    {
        log<INFO>("*** Rescaling solution with _scales");
        _scale_system(_Sxx,_x,_scales);
    }

    //double cond = N.cond();
    //   log<SAVE>("COND: ") % cond;

#ifdef DEBUG_SOLVER
    cerr<<"--- void Ls_neq::solve( type, bool ) "
            <<": "<<tim.toc()<<" s "<<endl;
#endif
}

// ...........................................................................
void Ls_neq::set_constraints(const ivg::Matrix & B_new,const ivg::Matrix & wgt_new,
                             const ivg::Matrix r_new)
// ...........................................................................
{
#ifdef DEBUG_SOLVER
    cerr<<"+++ void Ls_neq::set_constraints(const ivg::Matrix &, const ivg::Matrix &  ) "<<endl;
    tictoc tim;
    tim.tic();
#endif
    if (r_new.rows()==0)
    {
        _rB.resize(_N.rows(),1);
        _rB.zero();
    }
    else
        _rB = r_new;

    _B = B_new;
    _wB = wgt_new;

    _dN = _B.transpose()*_wB.diag()*_B;
    _dn = _B.transpose()*_wB.diag()*_rB;

#ifdef DEBUG_SOLVER
    cerr<<"--- void Ls_neq::set_constraints(const ivg::Matrix &, const ivg::Matrix &  )"
            <<": "<<tim.toc()<<" s "<<endl;
#endif
}

// ...........................................................................
void Ls_neq::update_constraints(const ivg::Matrix & B_new,const ivg::Matrix & wgt_new,
                                const ivg::Matrix r_new)
// ...........................................................................
{
#ifdef DEBUG_SOLVER
    cerr<<"+++ void Ls_neq::update_constraints(const ivg::Matrix &, const ivg::Matrix &  ) "<<endl;
    tictoc tim;
    tim.tic();
#endif
    ivg::Matrix rB = r_new;
    if (rB.rows()==0)
    {
        rB.resize(B_new.rows(),1);
        rB.zero();
    }

    // no pseudo observations until now
    if (_B.rows()==0)
    {
        _B = B_new;
        _wB = wgt_new;
        _rB = rB;

        _dN = B_new.transpose()*wgt_new.diag()*B_new;
        _dn = B_new.transpose()*wgt_new.diag()*rB;
    }
        // append new constraints to constraint matrices
    else
    {
        _B.append_rows(B_new);
        _wB.append_rows(wgt_new);
        _rB.append_rows(rB);

        _dN += B_new.transpose()*wgt_new.diag()*B_new;
        _dn += B_new.transpose()*wgt_new.diag()*rB;
    }


#ifdef DEBUG_SOLVER
    cerr<<"--- void Ls_neq::update_constraints(const ivg::Matrix &, const ivg::Matrix &  )"
            <<": "<<tim.toc()<<" s "<<endl;
#endif
}
// ...........................................................................
void Ls_neq::reduce_params(std::vector<int> idx)
// ...........................................................................
{
#ifdef DEBUG_SOLVER
    cerr<<"+++ void Ls_neq::reduce_param( std::vector<int> ) "<<endl;
    tictoc tim;
    tim.tic();
#endif
  
    _check_negative_index(idx);
    ;

    if (_scales.rows()!=0)
    {
        _scales.rem_r(idx);
    }

    // indexes of params which will not be reduced
    std::vector<int> idx0;
    vector<int>::iterator it;
    for (int i = 0; i<_N.rows(); ++i)
    {
        it = find(idx.begin(),idx.end(),i);
        if (it==idx.end())
            idx0.push_back(i);
    }

    // divide x in two parts: x1 remines, x2 is reduced
    // => divide N in 4 sub-matrices and n in two vectors

    
    ivg::Matrix N11 = _N.get_sub(idx0,idx0);
   
    //for(int i=0;i<idx.size();i++)
    //  std::cout << idx.at(i) <<" ";
    //std::cout << endl;
    ivg::Matrix N12 = _N.get_sub(idx0,idx);
    
    ivg::Matrix T = _N.get_sub(idx,idx0); // N21
    
    ivg::Matrix N22 = _N.get_sub(idx,idx);
    
    ivg::Matrix n1 = _n.get_sub(idx0,{0});
   
    ivg::Matrix n2 = _n.get_sub(idx,{0});


    
    if (_dN.rows()!=0)
    {
        // get constraints for reduced and kept parameters
        ivg::Matrix dN22 = _dN.get_sub(idx,idx);
        ivg::Matrix dn2 = _dn.get_sub(idx,{0});

        // entries between these groups should be zero
        ivg::Matrix dN12 = _dN.get_sub(idx0,idx);

        if (((dN12.norm()).norm())(0)!=0.0)
        {
            stringstream errormessage;
            errormessage<<"void Ls_neq::reduce_params(...): Fatal Error - off-digonal terms in dN"<<endl;
            throw runtime_error(errormessage.str());
        }

        // resizing matrices
        _dN = _dN.get_sub(idx0,idx0);
        _dn = _dn.get_sub(idx0,{0});

        // get indixes of all rows
        vector<int> all_idx(_B.rows());
        iota(begin(all_idx),end(all_idx),0); // fill from 0 to #_B.rows-1
        _B = _B.get_sub(all_idx,idx0);

        N22 += dN22;
        n2 += dn2;
    }

   
    ivg::Matrix t = n2;
    try
    {
        T.solve_neq(N22); // N22^-1 * N21
        t.solve_neq(N22); // N22^-1 * n2   
    }
    catch (std::runtime_error& e)
    {
        string error_message = e.what();
        throw runtime_error("void Ls_neq::reduce_params( std::vector<int> ): Failed to reduce parameters: "+error_message);
    }

    // reduce params
    _N = N11-N12*T;
    _n = n1-N12*t;
   
    // calculate new btpb
    ivg::Matrix tmp = n2.transpose()*t;
    _btPb = _btPb-tmp(0);

#ifdef DEBUG_SOLVER
    cerr<<"--- void Ls_neq::reduce_param( std::vector<int> )"
            <<": "<<tim.toc()<<" s "<<endl;
#endif
}
//...........................................................................
void Ls_neq::calc_correct_btPb()
//...........................................................................
{
#ifdef DEBUG_SOLVER
    cerr<<"+++ void Ls_neq::calc_correct_btPb()"<<endl;
    tictoc tim;
    tim.tic();
#endif
    // IN CASE OF DGF NOT CORRECT RIGHT NOW !!!! AND IAA !!!!
    // we need to know:  N  n  nobs  nparam  sig02
    double vtpv_new = _sigma0_post*(_nobs-_nparam);
    ivg::Matrix xT = _x.transpose();
    ivg::Matrix part1 = xT*_N*_x;
    ivg::Matrix part2 = (xT*2)*_n;
    double ltpl_corr = vtpv_new-part1(0)+part2(0);

    stringstream ss;
    ss<<setprecision(8);
    ss<<"*** NoO+NoC: "<<_nobs<<" / numPar: "<<_nparam<<" / VarianceFactor: "<<_sigma0_post<<endl;
    ss<<"*** vtpv_old: "<<_rtPr<<" / vtpv_new: "<<vtpv_new<<" / ltpl_old: "<<_btPb<<" / ltpl_new: "<<ltpl_corr;

    log<DETAIL>(ss.str());

    _btPb = ltpl_corr;

#ifdef DEBUG_SOLVER
    cerr<<"--- void Ls_neq::calc_correct_btPb()"
            <<": "<<tim.toc()<<" s "<<endl;
#endif
}
//...........................................................................
void Ls_neq::scale_system()
//...........................................................................
{
#ifdef DEBUG_SOLVER
    cerr<<"+++ Ls_neq::scale_system( vector<double> )"<<endl;
    tictoc tim;
    tim.tic();
#endif

    if (_scales.rows()!=0)
        _scale_system(_N,_n,_scales);
    else
        throw runtime_error("void Ls_neq::scale_system(vector<double>): Fatal Error - _scales not set. No scaling possible.");

#ifdef DEBUG_SOLVER
    cerr<<"--- void Ls_neq::scale_system( vector<double> )"
            <<": "<<tim.toc()<<" s "<<endl;
#endif  
}
//...........................................................................
void Ls_neq::_scale_system(ivg::Matrix &M,ivg::Matrix &m,ivg::Matrix &scales)
//...........................................................................
{
#ifdef DEBUG_SOLVER
    cerr<<"+++ Ls_neq::_scale_system(ivg::Matrix, ivg::Matrix, vector<double> )"<<endl;
    tictoc tim;
    tim.tic();
#endif
//       log<INFO>("*** _scale_system()");
//       scales.show();

    if (M.rows()!=scales.rows()||M.cols()!=scales.rows())
        throw runtime_error("void Ls_neq::_scale_system(ivg::Matrix, ivg::Matrix, vector<double> ): Fatal Error - unequal length of scales and matrix");
    else
    {
        ivg::Matrix F = scales.diag();
        M = F*M * F;
        m = m.mult_elem(scales);
    }

#ifdef DEBUG_SOLVER
    cerr<<"--- Ls_neq::_scale_system(ivg::Matrix, ivg::Matrix, vector<double> )"
            <<": "<<tim.toc()<<" s "<<endl;
#endif  
}
//...........................................................................
void Ls_neq::enlarge_neq(int num)
//...........................................................................
{
#ifdef DEBUG_SOLVER
    cerr<<"+++ void Ls_neq::enlarge_neq( int )"<<endl;
    tictoc tim;
    tim.tic();
#endif

    if (num>0)
    {
        ivg::Matrix N(_N.rows()+num,_N.cols()+num,0.0);
        N.set_sub(0,0,_N);
        _N = N;

        ivg::Matrix n_add(num,1,0.0);
        _n.append_rows(n_add);

        ivg::Matrix scale_add(num,1,0.0);
        if (_scales.rows()!=0)
            _scales.append_rows(scale_add);

    }

#ifdef DEBUG_SOLVER
    cerr<<"--- void Ls_neq::enlarge_neq( int )"
            <<": "<<tim.toc()<<" s "<<endl;
#endif
}
//...........................................................................
void Ls_neq::stack_neq(ivg::Ls_neq &other,vector<int> positions)
//...........................................................................
{
#ifdef DEBUG_SOLVER
    cerr<<"+++ void Ls_neq::stack_neq(ivg::Ls_neq &, vector<int> )"<<endl;
    tictoc tim;
    tim.tic();
#endif

    //   _N.show();
    //   other._N.show();

    _check_negative_index(positions);

    if (other._N.rows()!=positions.size()||other._N.cols()!=positions.size())
    {
        stringstream errormessage;
        errormessage<<"void Ls_neq::stack_neq(...): Fatal Error - unequal length"<<endl;
        throw runtime_error(errormessage.str());
    }
    else
    {

        ivg::Matrix add_N(_N.rows(),_N.cols(),0.0);
        ivg::Matrix add_n(_N.rows(),1,0.0);
        ivg::Matrix add_scales(_N.rows(),1,0.0);

        for (int i_other = 0; i_other<positions.size(); i_other++)
        {
            int i_this = positions.at(i_other);
            for (int j_other = 0; j_other<positions.size(); j_other++)
            {
                int j_this = positions.at(j_other);
                add_N(i_this,j_this) = other._N(i_other,j_other);
                add_N(j_this,i_this) = other._N(j_other,i_other);
            }

            add_n(i_this) = other._n(i_other);

            if (_scales.rows()!=0)
                _scales(i_this) = other._scales(i_other);
        }

        // final stacking
        _N += add_N;
        _n += add_n;
        _btPb += other._btPb;

    }


    //   abort();

#ifdef DEBUG_SOLVER
    cerr<<"--- void Ls_neq::stack_neq(ivg::Ls_neq &, vector<int> )"
            <<": "<<tim.toc()<<" s "<<endl;
#endif
}
//...........................................................................
void Ls_neq::epoch_transformation(int idx_offset,int idx_rate,double delta_t)
//...........................................................................
{
#ifdef DEBUG_SOLVER
    cerr<<"+++ void Ls_neq::epoch_transformation( int idx_offset, int idx_rate, double delta_t )"<<endl;
    tictoc tim;
    tim.tic();
#endif  

    ivg::Matrix TT;
    TT.eye(_nparam);

    TT(idx_offset,idx_rate) = delta_t;

    _N = TT*_N*TT.transpose();
    _n = TT * _n;

#ifdef DEBUG_SOLVER
    cerr<<"--- void Ls_neq::epoch_transformation( int idx_offset, int idx_rate, double delta_t )"
            <<": "<<tim.toc()<<" s "<<endl;
#endif
}
//...........................................................................
void Ls_neq::apriori_transformation(ivg::Matrix new_aprioris,ivg::Matrix old_aprioris)
//...........................................................................
{
#ifdef DEBUG_SOLVER
    cerr<<"+++ void Ls_neq::apriori_transformation( ivg::Matrix, ivg::Matrix )"<<endl;
    tictoc tim;
    tim.tic();
#endif  
   
    if(old_aprioris.rows() != new_aprioris.rows() || old_aprioris.cols() != new_aprioris.cols() || old_aprioris.cols()!=1 )
    {
        stringstream errormessage;
        errormessage << "void Ls_neq::apriori_transformation( ... ): Unequal Length of apriori vectors" << endl;
        throw runtime_error( errormessage.str() );
    }
    else
    {
        ivg::Matrix t = new_aprioris - old_aprioris;
        
        _n = _n - _N * t; 

        ivg::Matrix tT = t.transpose();
        ivg::Matrix tmp = tT * (_n*2 - _N * t);
        _btPb = _btPb - tmp(0);
    }
   
#ifdef DEBUG_SOLVER
    cerr<<"--- void Ls_neq::apriori_transformation(ivg::Matrix, ivg::Matrix )"
            <<": "<<tim.toc()<<" s "<<endl;
#endif
}
//...........................................................................
void Ls_neq::_check_negative_index(vector<int> &idx)
//...........................................................................
{
#ifdef DEBUG_SOLVER
    cerr<<"+++ void Ls_neq::_check_negative_index(vector<int> idx)"<<endl;
    tictoc tim;
    tim.tic();
#endif

    if (find(idx.begin(),idx.end(),-1)!=idx.end())
    {
        copy(idx.begin(),idx.end(),ostream_iterator <int>(cerr," | "));
        cerr<<endl;
        stringstream errormessage;
        errormessage<<"void Ls_neq::_check_negative_index(vector<int> idx): Negative Index -1 not possible"<<endl;
        throw runtime_error(errormessage.str());
    }


#ifdef DEBUG_SOLVER
    cerr<<"--- void Ls_neq::_check_negative_index(vector<int> idx)"
            <<": "<<tim.toc()<<" s "<<endl;
#endif
}
// ...........................................................................
void Ls_neq::show(bool verbose)
// ...........................................................................
{
    cerr << "++++++++++++ Ls_neq.show() +++++++++++++++++" << endl;
//    cout << " _start_epoch: " << setw(8) << left << setprecision(3) << fixed << _start_epoch.get_double_mjd() << endl;
    _N.cout_size(); cerr << "_N" << endl;
    if(verbose)
        _N.show();
    _n.cout_size(); cerr << "_n" << endl;
    if(verbose)
        _n.show();
    _dN.cout_size(); cerr << "_dN" << endl;
    if(verbose)
        _dN.show();
    
//    ivg::Matrix U,S,VT;
//    _N.svd( U, S, VT );
//    S.show();
    cerr << "------------ Ls_neq.show() -----------------" << endl;
    
}

void Ls_neq::get_neq(ivg::Matrix & N, ivg::Matrix & n, bool cnstr)
{
   N = _N;
   n = _n;
   
   if(cnstr)
   {
     if( _dN.rows() != _N.rows() || _dn.rows() != _n.rows()) {
       if (_dN.rows()>0 || _dn.rows()>0)
          throw runtime_error("void Ls_neq::get_neq(ivg::Matrix & N, ivg::Matrix & n, bool cnstr): Dimensions of _dN/_dn and _N/_n does not correspond.\n");
     }
      else
      {
        N += _dN; 
        n += _dn;
      }
   }
}

// ===========================================================================
// 		private methods
// ===========================================================================

//ivg::Matrix scaler(M.rows(),M.cols(),1.0);
//
//       //only scale if scale unequal 1.0 we need to scale
//       for(int i=0; i<scales.rows(); i++){
//
//           if(scales(i) != 1.0)
//           {
//               int row = i;
//               for(int c = 0; c < M.cols(); c ++ )
//               {
//                   scaler(row,c) = scales(i) * scales(i);
//               }
//
//               int column = i;
//               for(int r = 0; r< M.cols(); r ++ )
//               {
//                   scaler(r,column) = scales(i) * scales(i) ;
//               }
//
//               m(i) = m(i) / scales(i);
//
//           }
//
//        }
//
//       M = M / scaler;
//       scaler.show();

}
