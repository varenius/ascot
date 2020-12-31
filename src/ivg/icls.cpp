# include "icls.h"

namespace ivg
{

// ...........................................................................
Icls::Icls()					// default constructor
// ...........................................................................
{
  _alpha = 0.05;
}

// constructor with initial solution, design matrix level
// ...........................................................................
Icls::Icls(ivg::Matrix &A_, ivg::Matrix &l_, ivg::Matrix &invSigma_, ivg::Matrix &B_, ivg::Matrix &b_, ivg::Matrix &x_)		
// ...........................................................................
{
    // initialize some values
    _A = A_;
    _oc = l_;
    ivg::Matrix invSigma = invSigma_;
    _alpha = 0.05;			// set level of significance to 5 percent
    _B = B_;
    _b = b_;
    _x_icls = x_;

    // resize _N, _n and set all elements to zero
    _N.resize( _A.cols(), invSigma.rows() );	
    _n.resize( _A.cols(), 1 ) ;
    _N.zero(); 
    _n.zero();

    // from observation equations to normal equations
    _n = _A.transpose() * invSigma * _oc;								
    _N = _A.transpose() * invSigma * _A; 								

    // Covariance matrix of the parameters of the unconstrained system
    _Sigma_x_ols = _N;
    _Sigma_x_ols.inv() ;

    //solution of unconstrained system
    _x_ols = _Sigma_x_ols * _n;

    // LSI2QP: transform Least Squares problem to general Quadratic Program
    // H = 2 * ( _A' * inv_Sigma * _A + reg * eye(size(_A'*_A,1)));		_N = _A' * inv_Sigma * _A; 	reg = 0 ;
    // h = (-2 * _oc' * inv_Sigma * _A)';								_n = _oc' * inv_Sigma * _A;

    // set objective function matrix _C
    ivg::Matrix tmp = _N;
    tmp *= 2.0;
    _C = tmp;

    // set objective function vector _c
    tmp = _n;
    tmp *= -2.0;
    _c = tmp;

    // compute number of params m and constraints p
    calc_nconst();
    calc_nparam();
    
    _wrms = 0.0;
    _rms = 0.0;
    _vfac = 0.0;
}

// 2nd constructor with initial solution, design matrix level (for ICLS-MC)
// ...........................................................................
Icls::Icls(ivg::Matrix &N_, ivg::Matrix &A_, ivg::Matrix &l_, ivg::Matrix &invSigma_, ivg::Matrix &B_, ivg::Matrix &b_, ivg::Matrix &x_)		
// ...........................................................................
{
    // initialize some values
    _A = A_;
    _oc = l_;
    ivg::Matrix invSigma = invSigma_;
    _alpha = 0.05;			// set level of significance to 5 percent
    _B = B_;
    _b = b_;
    _x_icls = x_;

    // resize _N, _n and set all elements to zero
    _N.resize( _A.cols(), invSigma.rows() );	
    _n.resize( _A.cols(), 1 ) ;
    _N.zero(); 
    _n.zero();

    // normal equation system
    _N = N_ ;
    _n = _A.transpose() * _oc;								

    // Covariance matrix of the parameters of the unconstrained system
    _Sigma_x_ols = _N;
    _Sigma_x_ols.inv() ;

    //solution of unconstrained system
    _x_ols = _Sigma_x_ols * _n;

    // LSI2QP: transform Least Squares problem to general Quadratic Program
    // H = 2 * ( _A' * inv_Sigma * _A + reg * eye(size(_A'*_A,1)));		_N = _A' * inv_Sigma * _A; 	reg = 0 ;
    // h = (-2 * _oc' * inv_Sigma * _A)';								_n = _oc' * inv_Sigma * _A;

    // set objective function matrix _C
    ivg::Matrix tmp = _N;
    tmp *= 2.0;
    _C = tmp;

    // set objective function vector _c
    tmp = _n;
    tmp *= -2.0;
    _c = tmp;

    // compute number of params m and constraints p
    calc_nconst();
    calc_nparam();
    
    _wrms = 0.0;
    _rms = 0.0;
    _vfac = 0.0;    
}


// constructor with initial solution, normal equation level
// ...........................................................................
Icls::Icls(ivg::Matrix &N_, ivg::Matrix &n_, ivg::Matrix &B_, ivg::Matrix &b_, ivg::Matrix &x_)		
// ...........................................................................
{
    // initialize some values
    _alpha = 0.05;			// set level of significance to 5 percent
    _B = B_;
    _b = b_;
    _x_icls = x_;

    _n = n_ ;
    _N = N_ ;

    // Covariance matrix of the parameters of the unconstrained system
    _Sigma_x_ols = _N;
    _Sigma_x_ols.inv() ;

    //solution of unconstrained system
    _x_ols = _Sigma_x_ols * _n;

    // LSI2QP: transform Least Squares problem to general Quadratic Program
    // H = 2 * ( _A' * inv_Sigma * _A + reg * eye(size(_A'*_A,1)));		_N = _A' * inv_Sigma * _A; 	reg = 0 ;
    // h = (-2 * _oc' * inv_Sigma * _A)';								_n = _oc' * inv_Sigma * _A;

    // set objective function matrix _C
    ivg::Matrix tmp = _N;
    tmp *= 2;
    _C = tmp;

    // set objective function vector _c
    tmp = _n;
    tmp *= -2;
    _c = tmp;

    // compute number of params m and constraints p
    calc_nconst();
    calc_nparam();
    
    _wrms = 0.0;
    _rms = 0.0;
    _vfac = 0.0;    
}

// constructor using ivg::Lsa object
// ...........................................................................
Icls::Icls( ivg::Lsa * lsa, ivg::Matrix &B_, ivg::Matrix &b_, ivg::Matrix &x_ )		
// ...........................................................................
{ 
    _verbose = true;
    _alpha = 0.05;	// set level of significance to 5 percent
    _B = B_;
    _b = b_;
    _x_icls = x_;        
    _lsa = lsa;
    _lsa->get_neq( _N, _n, true );
    _A = *(_lsa->get_design_ptr());
    _oc = *(_lsa->get_obs_ptr());
    
    // Covariance matrix of the parameters of the unconstrained system
    _Sigma_x_ols = _N;
    _Sigma_x_ols.inv() ;

    //solution of unconstrained system
    _x_ols = _Sigma_x_ols * _n;

    // LSI2QP: transform Least Squares problem to general Quadratic Program
    // H = 2 * ( _A' * inv_Sigma * _A + reg * eye(size(_A'*_A,1)));		_N = _A' * inv_Sigma * _A; 	reg = 0 ;
    // h = (-2 * _oc' * inv_Sigma * _A)';								_n = _oc' * inv_Sigma * _A;
    
    // set objective function matrix _C
    ivg::Matrix tmp = _N;
    tmp *= 2.0;
    _C = tmp;

    // set objective function vector _c
    tmp = _n;
    tmp *= -2.0;
    _c = tmp;
    
    // compute number of params m and constraints p
    calc_nconst();
    calc_nparam();
    
    _wrms = 0.0;
    _rms = 0.0;
    _vfac = 0.0;    
}


// ==============================================
// =========== MEMBER-functions: ================
// ==============================================
// ...........................................................................
void Icls::set_statistics( double wrms, double rms, double vfac )
// ...........................................................................
{
    _wrms = wrms;
    _rms = rms;
    _vfac = vfac;
}

// ...........................................................................
void Icls::get_statistics( double & wrms, double & rms, double & vfac )
// ...........................................................................
{
    wrms = _wrms;
    rms = _rms;
    vfac = _vfac;
}

// ...........................................................................
ivg::Matrix Icls::get_confidence_limits()
// ...........................................................................
{
    ivg::Matrix conf_limits( _x_icls.size(1), 2 );
    conf_limits.set_sub( 0, 0, _x_icls - _cl_icls(":",0) );
    conf_limits.set_sub( 0, 1, _cl_icls(":",1) - _x_icls );
    
    return conf_limits;    
}

// compute number of constraints
// ...........................................................................
void Icls::calc_nconst()		
// ...........................................................................
{
  _p = _b.rows();
}

// compute number of parameters
// ...........................................................................
void Icls::calc_nparam()		
// ...........................................................................
{
  _m = _c.rows();
}

// solve quadratic problem with Active-Set-Method (according to Woelle)
// ...........................................................................
void Icls::solve_with_active_set()		
// ...........................................................................
{  
    if( _verbose )
        log<RESULT>("*** Start with Active Set ");
    
    int s;
    double qmin;
    ivg::Matrix bActive, W, w,  V, v, W_T, help, help2, W_, g, Proj, p_, Vjp, Vjz, q, Vneu, vneu, Wneu, wneu;
    ivg::Matrix H, Htmp, rhs, zeros, pk;

    // identify active constraints
    ivg::Matrix BT = _B.transpose() ;

    bActive = BT * _x_icls;
    bActive -= _b;


    // ************************************************************************************
    // find active and inactive constraints and subdivide _B in W (active) and V (inactive)
    // ************************************************************************************
    w.resize( 1, 1 );
    W.resize( _B.rows(), 1 );
    v.resize( 1, 1 );
    V.resize( _B.rows(), 1 );		

  
    
    ivg::Matrix wIdx(1,1,0.0);
    ivg::Matrix vIdx(1,1,0.0);

    for ( int i = 0; i < _p; i++)
    {
        if( bActive(i,0) == 0 )
        {			
            ivg::Matrix wNew(1,1);
            wNew = _b.get_sub( i,0,i,0 );
            w.append_rows( wNew );

            ivg::Matrix WNew( _B.rows(),1 );
            WNew = _B.get_sub( 0,i,_B.rows()-1,i );
            W.append_cols(WNew);					

            wIdx.append_rows( i );
        }
        else
        {
            ivg::Matrix vNew(1,1);
            vNew = _b.get_sub( i,0,i,0 );		
            v.append_rows( vNew );

            ivg::Matrix VNew( _B.rows(),1 );
            VNew = _B.get_sub( 0,i,_B.rows()-1,i );
            V.append_cols(VNew);		

            vIdx.append_rows( i );	
        }
    }

    w.rem_r( 0 );
    W.rem_c( 0 );
    v.rem_r( 0 );
    V.rem_c( 0 );
    wIdx.rem_r( 0 );
    vIdx.rem_r( 0 );
   

    // >>> Start of main loop <<<  
    // cf. algorithm 1: Active Set method, Roese-Koerner (2009)
    for ( int i = 0; i < 1000; i++)
    {
        // show current solution
//        if( _verbose )
//        {
//            cerr << i+1 << ". Iteration " << endl;
//            cerr << w.rows() << " active constraints." << endl; 
//            wIdx.show();
//            vIdx.cout_size();
//        }

        // check if no constraint is violated
        help2.resize(_B.cols(), 1);
        help2.zero();
        help2.plus_product_of( _B, _x_icls, CblasTrans, CblasNoTrans );

        for (int counter=0; counter < _b.rows(); counter++)
        {
            if (help2(counter,0) - _b(counter,0) > 1e-5)	
            {
                stringstream errormessage;
                errormessage << "Constraint is violated by "  << help2(counter,0)-_b(counter,0) << endl;
                throw runtime_error( errormessage.str() );
            }
        }

        // test maximal constraint violation
        if ( W.size(2) != 0 && w.size(1) != 0 )
        {	
            ivg::Matrix TestViolation = _B.transpose() * _x_icls - _b ;
            // cout << "Max constraint violation: " << TestViolation.max() << endl;
        }

        // *********************************************************************************
        // >>> calculate gradient g, projector, search direction p and increment q (scalar)
        // *********************************************************************************
 
        // gradient
        g= _C * _x_icls;					
        g = g + _c;
       
        if ( W.cols() > 0.0 )	// are there active constraints?			
        {
            // >>>> nach BEST (1984, p. 74 ) <<<<                     
            // 'classic' constrained estimation approach of the extended 
            // (normal) equation system
            H = _C ;
            H.append_cols( W );
            Htmp = W ;
            Htmp.transpose();
            zeros.resize( (int)(W.size(2)), (int)(W.size(2)), 0.0 ) ;
            Htmp.append_rows( zeros ) ;
            H.append_rows( Htmp.transpose() ) ;

            // right hand side
            rhs = g * -1.0 ;
            zeros.resize( (int)(W.size(2)), 1, 0.0 ) ;
            rhs.append_rows( zeros );

            // solve extended equation system
            pk = rhs ;
            pk.solve_geeqs( H ) ;
            p_ = pk.get_sub( 0, 0, (int)(W.size(1))-1, 0 ) ;
            _k1 = pk.get_sub( (int)(W.size(1)), 0, (int)(pk.size(1))-1, 0 ) ;          
        }
        else 	// if there are no active constraints
        {
            p_ = g * -1.0 ;
            p_.solve_neq( _C );	
        }
             
        if ( v.rows() > 1 )
        {
            // initialize with big number
            q.resize( V.rows(),1 );		
            ivg::Matrix qtmp( q.size(1),1,1e5 );
            q = qtmp;
        }
        else
        {
            // initialize with big number
            q.resize(1,1);
            q(0,0) = 1e5;
        }

        // find inactive constraints that might get violated
        for (int  j = 0; j < v.rows(); j++)
        {		
            ivg::Matrix VjSlice( V.rows(), 1 );
            VjSlice = V.get_sub( 0, j, V.rows()-1, j );

            ivg::Matrix Vj(VjSlice);
            Vjp.resize(1,1);					  	// set Vjp = [0]
            Vjp(0,0) = 0;		  
            Vjp.plus_product_of( Vj, p_, CblasTrans, CblasNoTrans );	//Vjp = Vj' *p_

            if ( Vjp(0,0) > 0 )
            {
                Vjz.resize(1,1);	// set Vjz = [0]
                Vjz(0,0) = 0;
                Vjz.plus_product_of( Vj, _x_icls, CblasTrans, CblasNoTrans ); // Vjz = Vj' *_x_icls

                // calculate distance to inactive constraint
                q(j,0) = (v(j,0) - Vjz(0,0)) / Vjp(0,0);	
            }
        }     
       
        // minimal value of q (increment; maximal step size)
        qmin = q.min(s);
         
//        if( _verbose )
//            cerr << "qmin (1)" << qmin << " at idx " << s << endl;         
         
         
        // ******************************************
        // distinction of cases (Fallunterscheidung) 
        // ******************************************
        if( qmin >= 1.0 )
        {
            // update solution vector
            _x_icls += p_;

            // check if no constraint is violated				
            if ( W.size(2) != 0 && w.size(1) != 0 )
            {	
                help2.resize(W.cols(),1);
                help2.zero();
                help2.plus_product_of( W, _x_icls, CblasTrans, CblasNoTrans ); 

                for ( int counter=0; counter < w.rows(); counter++ )
                {
                    if ( help2(counter, 0) - w(counter, 0) > 1e-5 )	
                    {
                        stringstream errormessage;
                        errormessage << "Former active constraint is violated by "  << help2(counter, 0)-w(counter, 0) << endl;
                        throw runtime_error( errormessage.str() );
                    }
                }
            }

            // update gradient
            g = _C * _x_icls + _c;			

            _k1.resize( _B.cols(), 1, 0.0 );					

            if ( (W.cols() > 0) )  		// are there former active constraints?
            {				
                if ( _k1.min() < 0.0 ) 
                {
                    for ( int j = 0; j < W.cols(); j++ )
                    {
                        if ( _k1(j, 0) < 0.0 ) 
                        {
                            if( _verbose )
                                log<RESULT>("*** Deactive active constraint");

                            // change constraint from W,w into V,v 										
                            ivg::Matrix vNew(1,1);
                            vNew = w.get_sub( j,0,j,0 );		
                            vNew.append_rows( v );
                            v = vNew;

                            ivg::Matrix VNew( W.rows(),1 );
                            VNew = W.get_sub( 0,j,W.rows()-1,j );
                            VNew.append_cols(V);
                            V = VNew ;

                            ivg::Matrix vIdxNew = wIdx.get_rows(j,j);
                            vIdxNew.append_rows( vIdx );
                            vIdx = vIdxNew ;
                            wIdx.rem_r(j);

                            W.rem_c(j);
                            w.rem_r(j);											
                        }
                    }
                }
                else
                {
                    // best possible solution has been successfully found
                    if( _verbose )
                        log<RESULT>("*** Solved with Active-Set-Method");

                    // set parameter vector in ls_solution
                    _x = _x_icls;
                    _Sxx = _Sigma_x_ols;
                    
                    // update vector of residuals in _lsa object 
                    // (needed for the calculation of RMS and WRMS)
//                      ivg::Matrix r = ( _oc - _A * _x_icls);                                                
//                    ivg::Matrix W = *(_lsa->get_wgt_ptr());
//                    ivg::Matrix Rinv = W.pow(-0.5);
//                    Rinv.to_diag_matrix();
//                    _r = Rinv * r;                            
                    
                    _r = ( _oc - _A * _x_icls); 
                    _lsa->set_resid( _r );
                    
                    _lambda.resize( _B.size(2),1,0.0 );
                    _lambda.setIdx( wIdx, _k1 );                    
                    _k1 = _k1.get_sub(0,0,w.rows()-1,0);
                    
                    _cnstr_idx = W.find_idx( ne, 0.0 );
                    
//                    if( _verbose )
//                        wIdx.show();
                
                    return;
                }
            }		 
        }
        else
        {
            if( _verbose )
                log<RESULT>("*** Insert new active constraint");

            // calculate new initial solution 
            p_ *= qmin;				// scale search direction
            _x_icls = _x_icls + p_;		// update parameter

            // change constraint linked with qmin from V,v into the active set W,w
            ivg::Matrix wNew(1,1);
            wNew = v.get_sub( s,0,s,0 );		
            wNew.append_rows( w );
            w = wNew ;

            ivg::Matrix WNew( V.rows(),1 );
            WNew = V.get_sub( 0,s,V.rows()-1,s );
            WNew.append_cols( W );
            W = WNew ;

            ivg::Matrix wIdxNew = vIdx.get_rows(s,s);
            wIdxNew.append_rows( wIdx );
            wIdx = wIdxNew ;
            vIdx.rem_r(s);

            V.rem_c(s);
            v.rem_r(s);		
        }			
    }	
}


// ...........................................................................
void Icls::estimate_quality( int M )
// ...........................................................................
{
    _verbose = false;
    
    // Monte Carlo for QP
    ivg::Matrix x_hat( _m, M, 0.0 );
    ivg::Matrix x_qp( _m, M, 0.0 );
    ivg::Matrix k_qp( _p, M, 0.0 );  
    ivg::Matrix n;
    
    // get cholesky decomposition of weight matrix
    ivg::Matrix R = _lsa->get_chol_wgtmat();
    
    // create random matrix (normal distribution)
    ivg::Matrix E( R.size(1), M, 0.0 );
    E.rand_norm( 0.0, 2.0e-11 );      
    ivg::Matrix dn = _A.transpose()* R * E;     
    
    // add some noise to the observation vector and create observations for M MC samples
    _lsa->get_neq( _N, _n, true );
    _dn.resize(0,0,0.0);
    _dn.repmat( _n, 1, M );
    _dn += dn;    
      
//    ivg::Matrix E( R.size(1), M, 0.0 );
//    E.rand_norm( 0.0, 2.0e-11 );     
//    ivg::Matrix dL = R * E;     
//    
//    // add some noise to the observation vector and create observations for M MC samples
//    _dL.repmat( *(_lsa->get_obs_ptr()), 1, M );
//    _dL += dL;  
        
    log<RESULT>( "*** Starting Monte-Carlo (" ) % M % " MC sweeps)" ;
    
    for ( int i=0; i < M; i++ )
    {
        if ( i % ( M / 10 ) == 0 )
            log<RESULT>("") % (i+1) % ". MC sweep";

        // update right-hand-side vectors n and c
        n = _dn( ":", i );
        update_MC_sample( n );
        
        solve_with_active_set();		
        x_qp.set_sub( 0, i, get_x() );
        k_qp.set_sub( 0, i, get_lagrange_mult() );				
    }

    // undo: set origninal o-c vector
    _lsa->get_neq( _N, _n, true );
    update_MC_sample( _n );
    _verbose = true;
    
    log<RESULT>( "*** Finished Monte-Carlo" );

    // compute histograms and confidence regions
    ivg::Matrix Sort_x_hat, Sort_x_qp;
    ivg::Matrix Conf_ols( x_qp.size(1), 2, 0.0 );   
    ivg::Matrix Conf_qp( x_qp.size(1), 2, 0.0 ); 

    // sort matrix by columns
    Sort_x_qp = x_qp;
    Sort_x_qp.sort_cols();
    Sort_x_hat = x_hat;
    Sort_x_hat.sort_cols();

    for ( int i = 0; i < x_qp.size(1); i++ )
    {
        Conf_ols(i,0) = Sort_x_hat( i, round( _alpha / 2 * M ) );			// lower bound OLS
        Conf_ols(i,1) = Sort_x_hat( i, round( ( 1 - ( _alpha / 2 ) ) * M )-1 );		// upper bound OLS
        Conf_qp(i,0) = Sort_x_qp( i, round( _alpha / 2 * M ) );				// lower bound ICLS
        Conf_qp(i,1) = Sort_x_qp( i, round( ( 1 - ( _alpha / 2 ) ) * M )-1 );		// upper bound ICLS
    }

    _cl_ols = Conf_ols;
    _cl_icls = Conf_qp;
    
}


// ...........................................................................
void Icls::update_MC_sample( ivg::Matrix &n )
// ...........................................................................
{
    // update right-hand-side of NEQ
//     ivg::Matrix n = _A.transpose() * l;

    // set objective function vector _c
    _c = n * -2.0;
     
    // update Jacobian matrix, o-c vector and initial parameter vector
    _A = *(_lsa->get_design_ptr());
    _oc = *(_lsa->get_obs_ptr());            
    ivg::Matrix one( _x_icls.rows(), 1, 1.0 );
    _x_icls = one;            

    // reset Lagrande multipliers
    _k1.resize( 1, 1, 0.0 );
    _lambda.resize( 1, 1, 0.0 );    
}     


} // namespace