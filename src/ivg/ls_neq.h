#ifndef LSNEQ_H
#define LSNEQ_H

#include "ls_solution.h"

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

namespace ivg {

    // transformation types 
    // cpwlf to offset or to offset+rate
    enum trafoto{ offset, offsetrate };

    // ===========================================================================

    class Ls_neq : public Ls_solution
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
         *  \return An instance of the class 'Lsa'
         */
        Ls_neq();

        /**
         *  \b Description: \n
         *  constructor using three input parameter: ivg::Matrix N, n, tx
         *  \param [in] [ivg::Matrix] Normal equation matrix N [n x n]
         *         [in] [ivg::Matrix] right hand side of Nx=n  [n x 1]
         *         [in] [ivg::Matrix] x-axis values of parameters x[1 x n]|[n x 1]: 
         *  \return An instance of the class 'Ls_neq'
         */
        Ls_neq(const ivg::Matrix& N, const ivg::Matrix& n, const int nparam, const ivg::Matrix& tx = ivg::Matrix());

        /**
         *  \b Description: \n
         *        Copy constructor
         *  \param [in] [ivg::Ls_neq] other observation
         *  \return An instance of the class 'Ls_neq'
         */
        Ls_neq(const ivg::Ls_neq &other);

        /**
         *  \b Description: \n
         *        Default deconstructor
         *  \param [in] no input parameters needed
         */
        ~Ls_neq();

        // ==============================================
        // =========== MEMBER-functions: ================
        // ==============================================

        /**
         *  \b Description: \n
         *  Method to get the number of observations
         *  \param [in] no input parameters needed
         *  \return [int] number of observations
         */
        int get_nobs() {
            return _nobs;
        };

        /**
         *  \b Description: \n
         *  Method to get the aposteriori variance factor
         *  \param [in] no input parameters needed
         *  \return [double] vfac
         */
        double calc_posterior_vfac() {
            return 999.0;
        };

        /**
         *  \b Description: \n
         *  Method to retrieve the post fit WRMS residuals: WRMS = sqrt((sum(r_i^2*p_i))/(sum p_i))
         *  where p_i are the weights according to the standard deviations of the observations
         *  \param [in] no input parameters needed
         *  \return [double] WRMS
         */
        double calc_wrms() {
            return 999.0;
        };

        /**
         *  \b Description: \n
         *  Method to retrieve the post fit RMS residuals: RMS = sqrt((sum(r_i^2*p_i))/(n))
         *  where is the number of observations 
         *  \param [in] no input parameters needed
         *  \return [double] RMS
         */
        double calc_rms() {
            return 999.0;
        };

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
        void trafo_params2polynomial(int deg, int idx, double t_new);

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
         *  Method to transform continuous picewise linear functions to something other.
         *  Other can be one (constant) parameter (type == offset) or to constant parameter + rate ( type == offsetrate )
         *  The cpwlf parameter will be removed and replaced by the new parameter.
         *  RIGHT NOW ONLY WORKING WITH n = 2 !!!
         *  \param [in] [ivg::trafoto] offset or offsetrate
         *         [in] [vector<int>] indexes of the CPWLF parameter with length n 
         *         [in] [vector<int>] time epochs of the CPWLF parameter with length n 
         *         [in] [double] t_0 epoch to be transformed to
         */
        void trafo_cpwlf2other(ivg::trafoto type, vector<int> idx, vector<double> t_i , double t_0);
        
        /**
         *  \b Description: \n
         *  Method to transform one (constant) parameter to additional linear parameter.
         *  Used for including stations velocities in global modus.
         *  \param [in] [vector<int>] column index according to constant basis parameter
         *         [in] [vector<double> ] time differences to be considered
         */
        void trafo_params2linear(vector<int> idx_0, vector<double> dt);
        
         /**
         *  \b Description: \n
         *  Method to detect unparameterized parameter in _neq_solution. This means zero in n_side as well as zeros
         *  in the related columns and rows of N. This is sometimes the case in e.g. nma2014a snx files.
         * 
         *  \return [vector<int>] which contains the origin indexes of the unparameterized parameters
         */
        vector<int> detect_erase_unparameterized();

        /**
         *  \b Description: \n
         *  Method to extend the matrix of constraint equations
         *  \param [in] [ivg::Matrix] linearized coeficients of r constraint equations [rxn]
         *         [in] [ivg::Matrix] vector of corresponding weights [rx1] 
         *         [in] [ivg::Matrix] right hand side of the constrained equations
         *                            (zero, if omitted)
         */
        void update_constraints(const ivg::Matrix & B_new,
                const ivg::Matrix & wgt_new,
                const ivg::Matrix r_new = ivg::Matrix());
        void update_nnt_nnr_constraints(const ivg::Matrix & B_new,
                const ivg::Matrix & wgt_new,
					 const ivg::Matrix r_new = ivg::Matrix(2, 2, 1.0))
         {
	   
	   update_constraints(B_new,wgt_new,r_new);
	 };
         void rm_nnt_nnr_constraints() {};
        /**
         *  \b Description: \n
         *  Method to set the matrix of constraint equations
         *  \param [in] [ivg::Matrix] linearized coeficients of r constraint equations [rxn]
         *         [in] [ivg::Matrix] vector of corresponding weights [rx1]
         *         [in] [ivg::Matrix] right hand side of the constrained equations
         *                            (zero, if omitted)
         */
        void set_constraints(const ivg::Matrix & B_new,
                const ivg::Matrix & wgt_new,
                const ivg::Matrix r_new = ivg::Matrix());

        /**
         *  \b Description: \n
         *  Method to build the normal equation system for a weighted least squares adjustment
         *  \param [in] [ivg::Matrix] functional/jacobian matrix A [m x n]
         *         [in] [ivg::Matrix] vector of observations (right hand side of Ax=b) [m x 1]
         *         [in] [ivg::Matrix] weight matrix/vector W [m x m]|[1 x m]|[m x 1]
         *         [in] [ivg::Matrix] x-axis values of parameters x[1 x n]|[n x 1]
         *         [in] [int] number of parameters if #param != n; e.g. jacobian matrix has been pre-reduced
         */
        void build_neq(const ivg::Matrix & A, const ivg::Matrix & W, const ivg::Matrix& b, const ivg::Matrix& epochs, const int nparam = 0);

        /**
         *  \b Description: \n
         *  Method to build the normal equation system for a least squares adjustment
         *  \param [in] [ivg::Matrix] functional/jacobian matrix A [m x n]
         *         [in] [ivg::Matrix] vector of observations (right hand side of Ax=b) [m x 1]
         *         [in] [ivg::Matrix] x-axis values of parameters x[1 x n]|[n x 1]
         *         [in] [int] number of parameters if #param != n; e.g. jacobian matrix has been pre-reduced
         */
        void build_neq(const ivg::Matrix & A, const ivg::Matrix& b, const ivg::Matrix& epochs, const int nparam = 0);

        /**
         *  \b Description: \n
         *  Method to set the elements of the normal equation system
         *  \param [in] [ivg::Matrix] Normal equation matrix N [n x n]
         *         [in] [ivg::Matrix] right hand side of Nx=n  [n x 1]
         *         [in] [ivg::Matrix] x-axis values of parameters x[1 x n]|[n x 1]: 
         *         [in] [int] number of parameters if #param != n; e.g. jacobian matrix has been pre-reduced
         */
        void set_neq(const ivg::Matrix & N, const ivg::Matrix & n, const ivg::Matrix& epochs, const int nparam = 0);

        /**
         *  \b Description: \n
         *  Method to get the normal equation matrix and vector either constrained
         *  or unconstrained
         *  \param [in] [bool] cnstr - return contrained (true) or unconstrained (false) system
         *  \return [ivg::Matrix] normal matrix
         *          [ivg::Matrix] normal vector
         */
        void get_neq(ivg::Matrix & N, ivg::Matrix & n, bool cnstr);

      void get_neq(ivg::Matrix & N, ivg::Matrix & n, ivg::Matrix & aplo, bool cnstr)
      {
	get_neq(N, n, cnstr);
	aplo=ivg::Matrix(n.rows(),n.cols(),0.0);
      }

        /**
         *  \b Description: \n
         *  Method to solve the normal equation system. 
         *  \param [in] [ivg::solutiontype] algorithm to be used
         *         [in] [bool] apply pre-conditioning prior to the solution (i.e. scale main diagonal of N to 1 and re-scale after the solution)
         */
        void solve(ivg::solutiontype type = ivg::solutiontype::neq_chol, bool pre_cond = true, double lambda=0.0);

        void new_solve();

        void reduce_params(std::vector<int> idx);

        void calc_correct_btPb();

        void enlarge_neq(int num);

        void scale_system();

        void set_scales(vector<double> scales) {
            _scales = ivg::Matrix(scales);
        };
        
        void clear_scales() { _scales.resize(0,0); };

        void stack_neq(ivg::Ls_neq &other, vector<int> positions);

        void epoch_transformation(int idx_offset, int idx_rate, double delta_t);

        void apriori_transformation(ivg::Matrix old_aprioris, ivg::Matrix new_aprioris);

        void _check_negative_index(vector<int> &idx);

        void show(bool verbose = false);

    private:
        // ==============================================
        // ============ private methods: ================
        // ==============================================
        void _scale_system(ivg::Matrix &M, ivg::Matrix &m, ivg::Matrix &scales);

        // ==============================================
        // ======== class variables / attributes: =======
        // ==============================================
        ivg::Matrix _N;
        ivg::Matrix _n;
        ivg::Matrix _dN;
        ivg::Matrix _dn;
        ivg::Matrix _scales;



    }; // ivg::Ls_neq
} // namespace 

#endif // LSNEQ_H
