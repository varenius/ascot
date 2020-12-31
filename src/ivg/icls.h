#ifndef ICLS_H
#define ICLS_H

#include "ls_solution.h"
#include "lsa.h"
#include "matrix.h"

/**
 *
 * @brief abstract class Icls 
 *         Inequality constrained least squares (ICLS) estimation deals with problems of the form
 *         
 *             min. (Ax - l)' (Ax - l)
 *             s.t. B'x <= b
 * 
 * @author SH - ivg::ASCOT developer team
 * @date 2016-01-19
 * @version 0.2
 *
 */

using namespace std;

namespace ivg {

// ===========================================================================
class Icls : public Ls_solution
// ===========================================================================
{

    public:
        // ==============================================
        // =============== Constructors: ================
        // ==============================================
        
        // default constructor
        Icls();					      
        
        // constructor with initial solution
        Icls(ivg::Matrix &A_, ivg::Matrix &l_, ivg::Matrix &invSigma_, ivg::Matrix &B_, ivg::Matrix &b_, ivg::Matrix &x_);
        
        // 2nd constructor with initial solution (for ICLS-MC)
	Icls(ivg::Matrix &N_, ivg::Matrix &A_, ivg::Matrix &l_, ivg::Matrix &invSigma_, ivg::Matrix &B_, ivg::Matrix &b_, ivg::Matrix &x_);
        
        // constructor using normal equations
        Icls(ivg::Matrix &N_, ivg::Matrix &n_, ivg::Matrix &B_, ivg::Matrix &b_, ivg::Matrix &x_);     
      
        // constructor using ivg::Lsa object
        Icls( ivg::Lsa * lsa, ivg::Matrix &B_, ivg::Matrix &b_, ivg::Matrix &x_ );             
        

        // ==============================================
        // =========== MEMBER-functions: ================
        // ==============================================
        
        /**
        *  \b Description: \n
        *        Method to update Monte Carlo sample: modify right-hand-side
        *        of the normal equation system (and the corresponding
        *        right-hand-side of the ICLS objective function matrix) 
        *  \param [in] new Monte Carlo sample: right-hand-side of NEQ system
        */        
        void update_MC_sample( ivg::Matrix &n );   

        /**
        *  \b Description: \n
        *        Method to update Monte Carlo sample: modify observation vector
        *  \param [in] new Monte Carlo sample: observation vector
        */        
        void setl( ivg::Matrix &l_ ){ _oc = l_; }
        
        void set_statistics( double wrms, double rms, double vfac );
        void get_statistics( double & wrms, double & rms, double & vfac );
        
        // from quadratic program
	ivg::Matrix get_x() const{ return( _x_icls ); };
	ivg::Matrix get_x_ols() const{ return( _x_ols ); };
	ivg::Matrix get_lagrange_mult() const{ return( _lambda ); };
        ivg::Matrix get_confidence_limits();
        std::vector<int> get_constr_indices(){ return _cnstr_idx; };
        
        ivg::Lsa * get_lsa_ptr(){ return _lsa; };
        	
	void calc_nconst();			// compute number of constraints p
	void calc_nparam();			// compute number of parameters m
	void solve_with_active_set();		// solve quadratic problem with Active-Set-Method (Woelle)
	void estimate_quality( int M );         // quality description via Monte Carlo simulation (M MC sweeps)        
        

    private:
        
        // ==============================================
        // ======== class variables / attributes: =======
        // ==============================================        
        // GMM, incl. ivg::Lsa object
        ivg::Matrix _A;			// design matrix
        ivg::Matrix _oc;		// observation vector
        ivg::Matrix _N;			// normal equation matrix
        ivg::Matrix _n;			// rhs of normal equations
        ivg::Matrix _Sigma;		// variance/covariance matrix of observations
        ivg::Matrix _Sigma_x_ols;	// variance/covariance matrix of parameters of ORDINARY LEAST SQUARES ESTIMATE (OLS)
        ivg::Matrix _r;	       
        ivg::Lsa * _lsa;                // ivg::Lsa object 
      
        // from quadratic program
	ivg::Matrix _C;			// objective function matrix
	ivg::Matrix _c;			// objective function vector
	ivg::Matrix _B;			// matrix of inequality constraints (transpose!!)
	ivg::Matrix _b;			// rhs of inequality constraints
	//ivg::Matrix _B_eq;		// matrix of equality constraints (transpose!!)
	//ivg::Matrix _b_eq;		// rhs of equality constraints
	ivg::Matrix _x_icls;		// vector of (estimated) parameters (constrained)
	ivg::Matrix _x_ols;		// vector of OLS solution
	ivg::Matrix _k1;		// vector of Lagrange Multipliers associated with inequalities
	ivg::Matrix _lambda;
	int _m;				// # parameters
	int _p;				// # inequalities
	//int _p_eq;			// # equalities
	//int _pw;			// # active inequalities  
        std::vector<int> _cnstr_idx;     // indices where inequality constraints are active
     
        // quantities needed for Monte-Carlo-VCV-estimatation 
        ivg::Matrix _X_icls;		// matrix containing the ICLS solutions x^(k) of different Monte-Carlos (MC) estimation steps
        ivg::Matrix _X_ols;             // matrix containing the OLS solutions x^(k) of different Monte-Carlos (MC) estimation steps
        ivg::Matrix _cl_icls;		// [m x 2] matrix containing confidence limits of the ICLS estimate
        ivg::Matrix _cl_ols;            // [m x 2] matrix containing confidence limits of the OLS estimate
        ivg::Matrix _dn;		// different realizations of the deviations of the observations
        double _alpha;                  // confidence limit              
        bool _verbose;                  // verbose flag to reduce the output in case of MC simulation
        
        double _wrms;
        double _rms;
        double _vfac;        
};

} // namespace 
#endif // ICLS_H
