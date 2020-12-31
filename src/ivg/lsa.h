#ifndef LSA_H
#define LSA_H

#include "ls_solution.h"
#include "ls_neq.h"
#include "math.h"
#include "boost/math/distributions/chi_squared.hpp"
#include "boost/math/distributions/students_t.hpp"
#include "boost/math/distributions/normal.hpp"
#include "boost/math/distributions/fisher_f.hpp"
#include "boost/math/distributions/chi_squared.hpp"

/**
 *
 * @brief abstract class lsa - for manipulating and solving classical least squares problem:
 *                       A*x = b
 *
 * @author TA - bakkari developer team
 * @date 2015-04-21
 * @version 0.2
 *
 */

using namespace std;
using namespace boost::math;

using boost::math::normal; // typedef provides default type is double
typedef students_t_distribution<double> students_t;
using boost::math::students_t;

namespace ivg {


    // ===========================================================================

    class Lsa : public Ls_solution
    // ===========================================================================
    {
        
    friend class Lsc;
        
    public:
        // ==============================================
        // =============== Constructors: ================
        // ==============================================
        /**
         *  \b Description: \n
         *        Default constructor
         *  \param [in] no input parameters needed
         *  \return An instance of the class 'Lsa'
         */
        Lsa();

        /**
         *  \b Description: \n
         *  constructor using three input parameter: ivg::Matrix A, b, W
         *  \param [in] [ivg::Matrix] functional/jacobian matrix A [n x m]
         *         [in] [ivg::Matrix] vector of observations (right hand side of Ax=b) [n x 1]
         *         [in] [ivg::Matrix] weight matrix/vector W [n x n]|[1 x n]|[n x 1]
         *         [in] [ivg::Matrix] x-axis values of observations b [1 x n]|[n x 1]: are set to equidistant spacing of starting at zero if none are given
         *         [in] [ivg::Matrix] x-axis values of parameters x[1 x m]|[m x 1]: are set to average value of tb if none are given
         *  \return An instance of the class 'Lsa'
         */
        Lsa(const ivg::Matrix& A, const ivg::Matrix b, const ivg::Matrix W,
                const ivg::Matrix tb = ivg::Matrix(), const ivg::Matrix tx = ivg::Matrix());
        
        /**
         *  \b Description: \n
         *        Default deconstructor
         *  \param [in] no input parameters needed
         */
        ~Lsa();

        // ==============================================
        // =========== MEMBER-functions: ================
        // ==============================================

        /**
         *  \b Description: \n
         *  Method to get the degrees of freedom
         *  \param [in] no input parameters needed
         *  \return [int] degrees of freedom
         */
        int get_degrees_of_freedom() {
            return _Qvv.trace();
        };

        /**
         *  \b Description: \n
         *  Method to get the number of observations
         *  \param [in] no input parameters needed
         *  \return [int] number of observations
         */
        int get_nobs() {
            return _A.rows();
        };

        /**
         *  \b Description: \n
         *  Method to retrieve the post fit residuals
         *  \param [in] no input parameters needed
         *  \return [ivg::Matrix] vector of residuals
         */
        ivg::Matrix get_resid();

        /**
         *  \b Description: \n
         *  Method to set the post fit residuals
         *  \param [in] [ivg::Matrix] vector of residuals
         */        
        void set_resid( ivg::Matrix r ){ _r = r; };

        /**
         *  \b Description: \n
         *  Method to set the cofactor matrix (CFM) of the post fit residuals
         *  \param [in] [ivg::Matrix] cofactor matrix (CFM) of the post fit residuals
         */                
        void set_CFM_resid( ivg::Matrix Qvv ){ _Qvv = Qvv; };

        void get_resid(ivg::Matrix &v, ivg::Matrix& Qvv);

        void get_wgt_matrix(ivg::Matrix& W, ivg::Matrix& factors);

        ivg::Matrix get_partial_redundancies() {
            return _Qvv.diag();
        };
        
        ivg::Matrix get_chol_wgtmat(){ return _R; };
        
        ivg::Matrix * get_design_ptr(){return &_A;};
        ivg::Matrix * get_design0_ptr(){return &_A0;};
        ivg::Matrix * get_obs_ptr(){return &_b;};
        ivg::Matrix * get_obs0_ptr(){return &_b0;};        
        ivg::Matrix * get_wgt_ptr(){return &_W;};
        ivg::Matrix * get_tobs_ptr(){return &_tb;};
        ivg::Matrix * get_aplo_ptr(){return &_aplo_vec;};       
        /**
         *  \b Description: \n
         *   Method to initialize an empty lsa object with dimension [n x m].
         *   Caution: only zero-matrices will be cretaed for _A, _b, _W and _tb, 
         *            _tx after filling the elements of these matrices, init-method
         *            has to be called to waranty a succefull estimation process!!!
         *  \param [in] [int] #observations
         *         [in] [int] #parameters
         *         [in] [bool] flag which indicates whether the covariance matrix will be 
         *                     full (true) or digonal (false; default)
         */
        void resize(int nobs,int nparam,bool full_vcm=false);

        /**
         *  \b Description: \n
         *   Method to be called after changing elements of the involved matrices
         *  \param [in] no input parameters needed
         */
        void reinit();
        
        /**
         *  \b Description: \n
         *  Method to get the adjusted observations
         *  \param [in] no input parameters needed
         *  \return [ivg::Matrix] adjusted observations
         */
        //      ivg::Matrix get_adjusted_obs(){ return _b0+_r; };

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
         *  where is the number of observations 
         *  \param [in] no input parameters needed
         *  \return [double] RMS
         */
        double calc_rms();

        /**
         *  \b Description: \n
         *  Method to get the aposteriori variance factor
         *  \param [in] no input parameters needed
         *  \return [double] vfac
         */
        double calc_posterior_vfac();

        /**
         *  \b Description: \n
         *  Method to eliminate one parameter
         *  \param [in] [int] column index according to parameter to be eliminated
         */
        void fix_param(std::vector<int> & idx);

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
                                      ivg::Matrix & aprioris );

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
			      const ivg::Matrix r_new = ivg::Matrix(2, 2, 1.0)); // r_new is always [0,0]!?
        void rm_nnt_nnr_constraints();
        /**
         *  \b Description: \n
         *  Method to export Jacobian matrix, observations and weights for use in solvebackend (just for testing)
         *  \param [in] [std::string] directory for export
         *         [in] [std::string] name for exported files
         *         [in] [std::vector<int>] indexes of parameters which should be exported (caution: all other parameters are assumed to be fixed)
         */
        void write_backend_gmm(std::string dir, std::string name,
                std::vector<int> param_idx);

        /**
         *  \b Description: \n
         *  Method to remove observations and parameters from A, b, W and B. Caution: normal equations, parameters, residuals and so on keep unchanged
         *  \param [in] [std::vector<int>] indexes of observations (rows of A) to be removed
         *         [in] [std::vector<int>] indexes of parameters (columns of A) to be removed
         */
        void remove_data(std::vector<int> rows, std::vector<int> cols);

        /**
         *  \b Description: \n
         *  Method observations from A, b, W.
         *  \param [in] [std::vector<int>] indexes of observations (rows of A) to be removed
         */
        void remove_observations( std::vector<int> rows );
        
        /**
         *  \b Description: \n
         *  Method to calculate condition of least squares problem (see Demmel (1997), Ch. 3.3)
         *  \param [in] no input parameters needed
         *  \return [double] condition
         */
        double cond_ls();

        /**
         *  \b Description: \n
         *  Method to solve the normal equation system. Depending on input parameters,
         *  the procedure is applied iteratively. In iteration #1 only gross errors are eliminated in further iterations quantile*sigma_i is used as threshold
         *  \param [in] [int] outlier_iterations (0=none, 1=only gross errors, >1 quantile*sigma_i & gross errors) (default: 2)
         *         [in] [double] threshold: thereshold for detecting gross errors [s] (default = 1ns)
         *         [in] [doudle] quantile of normal distribution for scaling standard deviations of residuals to define outlier threshold (default: 3.0)
         *         [in] [bool] apply pre-conditioning prior to the solution (i.e. scale main diagonal of N to 1 and re-scale after the solution)
         *         [in] [double] regularization parameter for Tikhonov regularization
         */
        void solve_neq(int outlier_iterations = 2, double threshold = 1e-9,
                       double quantile = 3.0, bool pre_cond = true, 
                       double lambda = 0.0);

        //void solve(int outlier_iterations = 2, double threshold = 1e-9,
        //        double quantile = 3.0, bool pre_cond = true);

        void reduce_params(std::vector<int> idx);

        void insert_break(int idx, double epoch);

        void combine_params(vector<int> idx, bool remove);

        void save_A_oc(string path) {
            _A.save_bin(path + "_A.bin");
            _b.save_bin(path + "_b.bin");
        }

        /**
         *  \b Description: \n
         *  Method to get the normal equation matrix and vector either constrained
         *  or unconstrained
         *  \param [in] [bool] cnstr - return contrained (true) or unconstrained (false) system
         *  \return [ivg::Matrix] normal matrix
         *          [ivg::Matrix] normal vector
         */
       ivg::Ls_neq get_neq(){ return _neq; };

        // get and set weight matrix/vector (needed to eliminate data before creating full populated weight matrices) 
        ivg::Matrix get_weights();
        void set_weights(ivg::Matrix wgt);


        /**
         *  \b Description: \n
         *  Method to detect outliers
         *  \param [in] [double] significance level
         *  \param [out] [double] percentage of detected outlier
         */
        void data_snooping(double alpha, std::string type, int & noutl, double & perc_outl, double & expect_perc_outliers);

        /**
         *  \b Description: \n
         *  Method to restore outliers (simple method)
         *  \param [in] [double] quantile
         *  \return [int] number of restored outliers
         */        
        int restore_outliers( double q );


        void get_neq(ivg::Matrix & N,ivg::Matrix & n,ivg::Matrix & aplo, bool cnstr);
        void get_neq(ivg::Matrix & N,ivg::Matrix & n,bool cnstr);

        /**
         *  \b Description: \n
         *  Method to build the normal equation system
         *  \param [in] no input parameters needed
         */
        void build_neq();
                
        /**
         *  \b Description: \n
         *  Method to find undefined parameters, i.e., only zeros in columns of
         *  _A and _B
         *  \param [in] no input parameters needed
         *  \return [std::vector<int>] indexes
         */
        std::vector<int> find_undefined_param();
        
        void repalce_param(const ivg::Matrix& x, const ivg::Matrix& Sxx);
        
    private:
        // ==============================================
        // ============ private methods: ================
        // ==============================================

        /**
         *  \b Description: \n
         *  Method to detect outliers called by Lsa::solve
         *  \param [in] [int] outlier_iterations (0=none, 1=only gross errors, >1 quantile*sigma_i & gross errors) (default: 2)
         *         [in] [double] threshold: thereshold for detecting gross errors [s] (default = 1ns)
         *         [in] [double] quantile of normal distribution for scaling standard deviations of residuals to define outlier threshold (default: 3.0)
         */
        void _detect_outliers(int iterations, double threshold = 1e-9,
                double quantile = 3.0);


        /**
         *  \b Description: \n
         *  Method to include the square root of the diagonal weights into jacobian matrix
         *  and vector of observations. In case of a full VCM, this is needed for the first solutions
         *  where outliers should be determined. 
         *  
         *  \param [in] no input parameters needed
         */
        void _modify_gmm_diagonal();

        double _get_degrees_of_freedom();

        ivg::Matrix _get_chol_wgtmat_inv();

        std::vector<int> _get_idx_wo_outliers();
        
        // ==============================================
        // ==============================================
        ivg::Matrix _A;  // jacobian matrix (after full decorrelation)
        ivg::Matrix _b;  // vector of observations after full decorrelation
        ivg::Matrix _A0; // jacobian matrix
        ivg::Matrix _b0; // vector of observations
        ivg::Matrix _W; // weight matrix/vector
        ivg::Matrix _R; // cholesky decomposition of weight matrix

        ivg::Matrix _aplo_vec;
      
        ivg::Matrix _tb; // x-axis values (epochs) of observations b

        ivg::Matrix _r; // post fit residuals
        ivg::Matrix _Qvv;
        ivg::Matrix _dQvv; // update of Qvv if params have been reduced

        double _wrms; // WRMS of postfit residuals
        double _rms; // RMS of postfit residuals

        bool _full_vcm; // flag which is true if VCM of obs = W^-1 has off-diagonal entries

        ivg::Ls_neq _neq;

        std::vector<int> _idx_outliers; // indexes of outliers
        ivg::Matrix _A_outliers; // rows of the jacobian matrix corresponding to outliers
        ivg::Matrix _b_outliers; // outliers

        std::vector<int> _rows_nnt_nnr; // rows for nnt/nrr constraints
    }; // ivg::Lsa
} // namespace 

#endif // LSA_H
