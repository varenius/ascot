#ifndef LSC_H
#define LSC_H

#include "ls_solution.h"
#include "ls_neq.h"
#include "math.h"
#include "lsa.h"

/**
 * @brief abstract class lsc - for least squares collocation applications:
 *
 * @author SH - ivg::ASCOT developer team
 * @date 2016-05-05
 * @version 0.2
 *
 */

using namespace std;

namespace ivg {

    enum stoch_param{ tropo, clock };

//    class Lsa; // forward declaration
    
    // ===========================================================================

    class Lsc : public Ls_solution
    // ===========================================================================
    {
        
    friend class Lsa;
        
    public:
        // ==============================================
        // =============== Constructors: ================
        // ==============================================
        /**
         *  \b Description: \n
         *        Default constructor
         *  \param [in] no input parameters needed
         *  \return An instance of the class 'Lsc'
         */
        Lsc();

        /**
         *  \b Description: \n
         *  constructor using a ivg::Lsa object
         *  \param [in] [ivg::Lsa *] lsa object 
         *  \return An instance of the class 'Lsc'
         */
        Lsc( ivg::Lsa * lsa, ivg::Matrix B = ivg::Matrix() );
        
        /**
         *  \b Description: \n
         *        Default deconstructor
         *  \param [in] no input parameters needed
         */
        ~Lsc();

        
        // ==============================================
        // =========== MEMBER-functions: ================
        // ==============================================

        /**
         *  \b Description: \n
         *  Method to get the pointer to the ivg::Lsa object
         *  \param [in] no input parameters needed
         *  \return pointer to the ivg::Lsa object
         */                
        ivg::Lsa * get_lsa_ptr(){ return _lsa; };
        
        /**
         *  \b Description: \n
         *  Method to get the number of observations
         *  \param [in] no input parameters needed
         *  \return [int] number of observations
         */
        int get_nobs() {
            return _lsa->get_design_ptr()->rows();
        };        

        /**
         *  \b Description: \n
         *  Method to get the number of parameters (deterministic and stochastic)
         *  \param [in] no input parameters needed
         *  \return [int] number of parameters (deterministic and stochastic)
         */
        int get_nparams() { return ( _x.rows() + _y.rows() ); };               
        
        /**
         *  \b Description: \n
         *  Method to get the degrees of freedom
         *  \param [in] no input parameters needed
         *  \return [int] degrees of freedom
         */
        int get_degrees_of_freedom() { return _lsa->get_degrees_of_freedom(); };  
        
        /**
         *  \b Description: \n
         *  Method to get the "assignment map" G which consist of a 
         *  matrix Gi for each station
         *  of the stochastic parameters
         *  \param [in] no input parameters needed
         *  \return [map<string,ivg::Matrix>] assignment map
         */        
        std::vector< ivg::Matrix > get_G(){ return _G; };        
        
        /**
        *  \b Description: \n
        *        Method to find the number of stochastic parameters per station
        *  \param [in] [int] number of stochastic parameters
        */
        void find_n_stoch_params();
        
        /**
        *  \b Description: \n
        *        Method to find the number of stochastic parameters per station
        *  \param [in] [int] number of stochastic parameters
        */
        std::vector<int> get_n_stoch_params(){ return _nsparam; };
        
        /**
         *  \b Description: \n
         *  Method to get the vector of stochastic parameter estimates
         *  \param [in] no input parameters needed
         *  \return [ivg::Matrix] stochastic parameter vector
         */
        ivg::Matrix get_stoch_params() { return _y; };

        /**
         *  \b Description: \n
         *  Method to get the variance covariance matrix of the 
         *  stochastic parameter estimates
         *  \param [in] no input parameters needed
         *  \return [ivg::Matrix] VCM of stochastic parameters
         */
        ivg::Matrix get_stoch_vcm() { return _Syy; };        
        
        /**
         *  \b Description: \n
         *  Method to solve the least squares collocation method
         *  \param [in] no input parameters needed
         */        
        void solve( );
        
        /**
         *  \b Description: \n
         *  Method to calculate the design matrix including the partial 
         *  derivatives for the stochastic parameters
         *  \param [in] [std::vector<int>] vector of indices indicating which 
         *                                 columns of the Jacobian matrix have 
         *                                 to be deleted and added to the 
         *                                 stochastic Jacobian matrix
         */        
        void calc_stoch_design( std::vector<int> idx );

        /**
         *  \b Description: \n
         *  Method to calculate the variance covariance matrix 
         *  of the stochastic parameters
         *  \param [in] no input parameters needed
         */        
        void calc_stoch_VCM( ivg::Matrix Qy );

        /**
         *  \b Description: \n
         *  Method to set the "assignment map" G which consist of a 
         *  matrix Gi for each station
         *  of the stochastic parameters
         *  \param [in] [map<string,ivg::Matrix>] assignment map
         */        
        void set_G( std::vector< ivg::Matrix > G );
        
        /**
         *  \b Description: \n
         *  Method to transform add new constraint equations.
         *  \param [in] [int] linearized constrained equations (i.e., update of the 
         *                    jacobian matrix)
         *         [in] [ivg::Matrix] corresponding weights 
         *         [in] [ivg::Matrix] right hand side of the constrained equations
         *                            (zero, if omitted)
         */
        void update_constraints(const ivg::Matrix & B_new,
                const ivg::Matrix & wgt_new,
                const ivg::Matrix r_new = ivg::Matrix(2, 2, 1.0)); // r_new is always [0,0]!?        

         void update_nnt_nnr_constraints(const ivg::Matrix & B_new,
                const ivg::Matrix & wgt_new,
					 const ivg::Matrix r_new = ivg::Matrix(2, 2, 1.0))
         {
	   
	   update_constraints(B_new,wgt_new,r_new);
	 };
         void rm_nnt_nnr_constraints() {};

	   // r_new is always [0,0]!? 
        
         /**
         *  \b Description: \n
         *  Method to reduce parameters according to 
         *  Niemeier (2002) 'Ausgleichungsrechnung', pp 286-289 and
         *  Funcke (1982) 'Verfahren zur Parameterelimination im GMM', AVN, 89, pp 112-122.
         *  \param [in] [std::vector<int>] vector of indices of parameters to be deleted
         */               
        void reduce_params( std::vector<int> idx );
        
        /**
         *  \b Description: \n
         *  Method to fix parameters
         *  \param [in] [int] column index according to the parameters to be eliminated
         */
        void fix_param( std::vector<int> & idx );
  
        /**
         *  \b Description: \n
         *  Method to handle stochastic parameters: create stochastic Jacobian
         *  matrix B and fix corresponding column in the original Jacobian matrix
         *  \param [in] [int] column index according to parameters 
         *                    which are changed to be stochastic
         */
        void handle_stoch_param( std::vector<int> & idx );
        
        /**
         *  \b Description: \n
         *  Method to retrieve the post fit WRMS residuals: WRMS = sqrt((sum(r_i^2*p_i))/(sum p_i))
         *  where p_i are the weights according to the standard deviations of the observations
         *  \param [in] no input parameters needed
         *  \return [double] WRMS
         */
        double calc_wrms();

        /**
         *  \b Description: \n
         *  Method to retrieve the post fit RMS residuals: RMS = sqrt((sum(r_i^2*p_i))/(n))
         *  where n is the number of observations 
         *  \param [in] no input parameters needed
         *  \return [double] RMS
         */
        double calc_rms();        

        /**
         *  \b Description: \n
         *  Method to get the a posteriori variance factor
         *  \param [in] no input parameters needed
         *  \return [double] a posteriori variance factor (vfac)
         */
        double calc_posterior_vfac();
        
        /**
         *  \b Description: \n
         *  Method to transform one (constant) parameter to higher degree polynomials.
         *  It is assumed that the basis parameter at position idx is of degree 0.
         *  Parameters of degree 1 to deg are inserted at the end of the parameter list.
         *  \param [in] [int] maximal degree of new polynomial coefficients
         *         [in] [int] column index according to parameter of degree 0
         *         [in] [double] new parameter epoch
         */
        void trafo_params2polynomial( int deg, int idx, double t_new, 
                                      ivg::Matrix & aprioris, ivg::Matrix & new_apriori );        
        
        /**
         *  \b Description: \n
         *  Method to transform one (constant) parameter to continuous picewise linear functions.
         *  The basis parameter will be removed and the new CPWLF parameters appear at the end of the parameter list.
         *  \param [in] [int] column index according to constant basis parameter
         *         [in] [ivg::Matrix] epochs for CPWLF parameters
         */
        void trafo_params2cpwlf(int idx, ivg::Matrix t_new);
    
        /**
         *  \b Description: \n
         *  Method to detect outliers
         *  \param [in] [double] significance level
         *  \param [out] [double] percentage of detected outlier
         */
        void data_snooping( double alpha, std::string type, int & noutl, double & perc_outl, double & expect_perc_outliers );        
        
    private:
        
        // ==============================================
        // ============ private methods: ================
        // ==============================================

             
        // ==============================================
        // ========== private attributes: ===============
        // ==============================================
        
        ivg::Lsa * _lsa; // ivg::Lsa object 
        ivg::Matrix _B;  // jacobian matrix of stochastic parameters (after full decorrelation)
        ivg::Matrix _B0; // jacobian matrix of stochastic parameters
        ivg::Matrix _Qy; // covariance matrix for stochastic parameters y
        ivg::Matrix _Qw; // covariance matrix for stochastic parameters w
        
//        ivg::Matrix _x;    // estimated deterministic parameters
//        ivg::Matrix _Sxx;  // covariance matrix of estimated deterministic parameters        
        ivg::Matrix _y;    // estimated stochastic parameters
        ivg::Matrix _Syy;  // covariance matrix of estimated stochastic parameters        
        
        ivg::Matrix _v;  // post fit residuals (after full decorrelation)
        ivg::Matrix _v0; // post fit residuals
        ivg::Matrix _Qv; // covariance matrix of post fit residuals (after full decorrelation)
        
        std::vector< ivg::Matrix > _G; // collocation assignment matrix G (per station))
        std::vector<int> _nsparam;     // number of stochastic parameters per station
        double _scale;  // scaling factor for the covariance matrix for stochastic parameters y
        
        // coefficients of the auto-correlation function (from Titov, 2000)
        double _alpha_clock = 2.64;
        double _beta_clock = 8.64;
        double _phi_clock = 0.33;
        double _alpha_tro = 6.24;
        double _beta_tro = 6.48;
        double _phi_tro = 0.82;        

    }; // ivg::Lsc
} // namespace 

#endif // LSC_H
