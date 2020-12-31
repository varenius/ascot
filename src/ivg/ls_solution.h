#ifndef LSSOL_H
#define LSSOL_H

#include "matrix.h"
#include "tictoc.h"
#include "logger.h"

/**
 *
 * @brief abstract class ls_solution - base class
 *
 * for manipulating and solving least squares problems:
 *                    the equation
 *                       A*x = b
 *                    is solved by minimizing
 *                       (A*x-b)^T*(A*x-b)
 *
 * @author TA - bakkari developer team
 * @date 2015-04-21
 * @version 0.2
 *
 */

using namespace std;

namespace ivg {

    // Solution Type

    enum solutiontype {
        qr, svd, neq_lu, neq_chol
    };

    // ===========================================================================

    class Ls_solution
    // ===========================================================================
    {
    public:
        // ==============================================
        // =============== Constructors: ================
        // ==============================================
        /**
         *  \b Description: \n
         *        Default constructor
         *  \param [in] no input parameters needed
         *  \return An instance of the class 'Ls_sol'
         */
        //      Ls_solution();

        /**
         *  \b Description: \n
         *        Default deconstructor
         *  \param [in] no input parameters needed
         */
        //      virtual  ~Ls_solution();


        // ==============================================
        // =========== MEMBER-functions: ================
        // ==============================================

        /**
         *  \b Description: \n
         *        Copy constructor
         *  \param [in] [ivg::Ls_solution] other observation
         *  \return An instance of the class 'Ls_solution'
         */
        //       Ls_solution( const ivg::Ls_solution &other ){};

        /**
         *  \b Description: \n
         *  Method to solve theequation system. Depending on input parameter, the least squares problem is solved
         *  directly (QR decomposition, or SVD) or via normal equations (LU or Cholosky decomposition).
         */
        //      virtual void solve_ls( ivg::solutiontype type ) = 0;

        /**
         *  \b Description: \n
         *  Method to get the aposteriori variance factor
         *  \param [in] no input parameters needed
         *  \return [double] vfac
         */
        virtual double calc_posterior_vfac() {
        };

        /**
         *  \b Description: \n
         *  Method to transform add new constraint equations.
         *  \param [in] [int] linearized constrained equations (i.e., update of the 
         *                    jacobian matrix)
         *         [in] [ivg::Matrix] corresponding weights 
         *         [in] [ivg::Matrix] right hand side of the constrained equations
         *                            (zero, if omitted)
         */

        virtual void update_nnt_nnr_constraints(const ivg::Matrix & B_new,
			      const ivg::Matrix & wgt_new,
						const ivg::Matrix r_new = ivg::Matrix()) {}; // r_new is always [0,0]!?
      
      virtual void rm_nnt_nnr_constraints() {};
      
        virtual void update_constraints(const ivg::Matrix & B_new,
                const ivg::Matrix & wgt_new,
					const ivg::Matrix r_new = ivg::Matrix()) {};

        /**
         *  \b Description: \n
         *  Method to eliminate one parameter
         *  \param [in] [int] column index according to parameter to be eliminated
         */
        virtual void fix_param(std::vector<int> & idx) {
        };

        /**
         *  \b Description: \n
         *  Method to handle deterministic and stochastic parameters
         *  (e.g., for filter estimation or least squares collocation)
         *  \param [in] [int] column index according to parameter to be eliminated
         */
        virtual void handle_stoch_param(std::vector<int> & idx) {
        };        
        
        /**
         *  \b Description: \n
         *  Method to reduce one parameter
         *  \param [in] [int] column index according to parameter to be reduced
         */
        virtual void reduce_params(std::vector<int> idx) {
        };

        /**
         *  \b Description: \n
         *  Method to calculate the wrms
         */
        virtual double calc_wrms() {
        };

        /**
         *  \b Description: \n
         *  Method to calculate the rms
         */
        virtual double calc_rms() {
        };

        /**
         *  \b Description: \n
         *  Method to transform one (constant) parameter to higher degree polynomials.
         *  It is assumed that the basis parameter at position idx is of degree 0.
         *  Parameters of degree 1 to deg are inserted at the end of the parameter list.
         *  \param [in] [int] maximal degree of new polynomial coefficients
         *         [in] [int] column index according to parameter of degree 0
         *         [in] [double] new parameter epoch
         */
        virtual void trafo_params2polynomial( int deg, int idx, double t_new,
                                              ivg::Matrix & aprioris ) {
        };

        virtual void trafo_params2polynomial( int deg, int idx, double t_new,
                                              ivg::Matrix & aprioris, 
                                              ivg::Matrix & new_apriori ) {
        };      
        
        /**
         *  \b Description: \n
         *  Method to transform one (constant) parameter to continuous picewise linear functions.
         *  The basis parameter will be removed and the new CPWLF parameters appear at the end of the parameter list.
         *  \param [in] [int] column index according to constant basis parameter
         *         [in] [ivg::Matrix] epochs for CPWLF parameters
         */
        virtual void trafo_params2cpwlf(int idx, ivg::Matrix t_new) {
        };

        /**
         *  \b Description: \n
         *  Method to get the number of parameter
         *  \param [in] no input parameters needed
         *  \return [int] number of parameter
         */
        int get_nparams() {
            return _nparam;
        };

        /**
         *  \b Description: \n
         *  Method to get the vector of parameter estimates
         *  \param [in] no input parameters needed
         *  \return [ivg::Matrix] parameter vector
         */
        ivg::Matrix get_parameters() {
            return _x;
        };

        /**
         *  \b Description: \n
         *  Method to get the variance covariance matrix of the parameter estimates
         *  \param [in] no input parameters needed
         *  \return [ivg::Matrix] VCM of parameters
         */
        ivg::Matrix get_vcm() {
            return _Sxx;
        };

        /**
         *  \b Description: \n
         *  Method to get the correlation matrix between the estimated parameters
         *  \param [in] no input parameters needed
         *  \return [ivg::Matrix] correlation matrix of parameters
         */
        ivg::Matrix get_correlation_matrix() {
            ivg::Matrix S = ((_Sxx.diag())^(-.5)).diag();
            return ( S * _Sxx * S);
        };

        virtual int get_degrees_of_freedom() {
        };

        virtual int get_nobs() {
        };

        double get_btPb() {
            return _btPb;
        };

        void set_vcm(ivg::Matrix &Sxx) {
            _Sxx = Sxx;
        };

        void set_rtPr(double rtPr) {
            _rtPr = rtPr;
        };

        void set_btPb(double btPb) {
            _btPb = btPb;
        };

        void set_nobs(int nobs) {
            _nobs = nobs;
        };

        void set_nparam(int nparam) {
            _nparam = nparam;
        };

        void set_sigma0_post(double vf) {
            _sigma0_post = vf;
        };

        void set_sigma0_pre(double vf) {
            _sigma0_pre = vf;
        };

        virtual void show() {
        };

        /**
         *  \b Description: \n
         *  Method to get the normal equation matrix and vector either constrained
         *  or unconstrained
         *  \param [in] [bool] cnstr - return contrained (true) or unconstrained (false) system
         *  \return [ivg::Matrix] normal matrix
         *          [ivg::Matrix] normal vector
         */
        virtual void get_neq(ivg::Matrix & N,ivg::Matrix & n,bool cnstr) {
        };
        virtual void get_neq(ivg::Matrix & N,ivg::Matrix & n,ivg::Matrix & aplo,bool cnstr) {
	};
    protected:
        // ==============================================
        // ============ private methods: ================
        // ==============================================

        // ==============================================
        // ======== class variables / attributes: =======
        // ==============================================
        ivg::Matrix _x; // estimated parameters
        ivg::Matrix _Sxx; // covariance matrix of estimated parameters
        ivg::Matrix _tx; // x-axis values (epochs) of parameters x

        int _nobs; // number of observations
        int _nparam; // number of original (possibly reduced) estimated parameters

        double _sigma0_pre; // a proiori variance factor
        double _sigma0_post; // a posteriori variance factor

        double _btPb; // squared sum of the observations
        double _rtPr; // squared sum of postfit residuals

        ivg::Matrix _B; // matrix of pseudo observations
        ivg::Matrix _rB; // right hand side of pseudo observations
        ivg::Matrix _wB; // and corresponding weights

    }; // class Ls_solution

} // namespace
#endif // LSSOL_H
