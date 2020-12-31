
/* 
 * File:   network.h
 * Author: corbin
 *
 * Created on 12. September 2017, 12:27
 */

#ifndef NETWORK_H
#define NETWORK_H

#include "triangle.h"
#include "matrix.h"
#include "session.h"
#include "baseline.h"
#include "logger.h"

#include <map>
#include <vector>
#include <cmath>

#ifdef USE_AHC
//for ahc
#include "dataanalysis.h"
#include "stdafx.h"
#endif

class Network {
public:
    
    /**
     *  \b Description: \n
     *        Default constructor
     *  \param [in] no input parameters needed
     *  \return An instance of the class 'Network'
     */
    Network();
    
    /**
     *  \b Description: \n
     *        constructor using pointer to a Session
     *  \param [in] [ivg::Session*] pointer to a Session
     *  \return An instance of the class 'Network'
     */
    Network(ivg::Session* session);

     /**
      *  \b Description: \n
      *        Default deconstructor
      *  \param [in] no input parameters needed
      */
    ~Network();

    enum Algorithem {
        manual, mn, ahc
    };

    /**
     *  \b Description: \n
     * For every baseline in the  _baselines map the integer ambiguities are computed
     * and stored in _baselines.shift
     * afterwars _shift_residuals_and_delays() is called
     * 3 Algorithems are available:
     *  (1) manual: usees the statistics ui to shift them manually
     *  (2) mn: most neighbours -> residual with most neighbours is reference
     *  (3) ahc: agglomerativ hirachical clustering -> mean of cluster with most residuals is reference
     * 
     *  \param [in] [Algorithem] 
     * 
     */
    void compute_baselinewise_integer_ambiguities(Algorithem alg, ivg::Session* session = nullptr);

     /**
     *  \b Description: \n
     *   if at least one triangle exists (min 3 baselines) the dependend baseline
     *   is shifted so that the loop closure is minimal
     * 
     *  \param [in] no input parameters needed
     * 
     */
    void apply_closure_condition();

     /**
     *  \b Description: \n
     * These three private attributes are initailized by this function:
     *   std::vector<short> _all_integerAmbiguities;
     *   ivg::Matrix _all_group_delays;
     *   ivg::Matrix _all_group_delays_sigmas;
     * 
     *  The baselinewise stored values are copied into one matrix/vector 
     *  The size of the matrix/vector correspondes to the size of the ncfiles
     *  They are sorted like the ncfile
     * 
     *  \param [in] no input parameters needed
     */
    void fill_all_delay_and_integer_matrix();

    
     /**
     *  \b Description: \n
     * computes the ionospheric correction
     * 
     *  \param [in] [Network] X-band, S-band
     *         [out] [ivg::Matrix] ionospheric correction, sigma of ionospheric correction
     *               [std::vector<short>] error_flag information (  0=OK , -1 = Missing, -2 = bad )
     */
    static void get_ionospheric_correction(const Network& X_net, const Network& S_net,
            ivg::Matrix& delta_tau_x, ivg::Matrix& delta_tau_x_sigma, std::vector<short>& error_flag);
    
    void update_resids_in_session(ivg::Session* session);

    /**
     *  \b Description: \n
     * prints Network related information to the console
     * 
     *  \param [in] no input parameters needed
     */
    void print_Network() const;
    
    void save_resids( const std::string file ) const;

  /**
     *  \b Description: \n
     * returns a vector containing the integer ambiguities (same order as in ncfile)
     * if no  integer ambiguity could be computed (e.g. missing observatioon due to useflag)
     * the value is 666 
     * make sure that fill_all_delay_and_integer_matrix is called befor using this getter 
     * 
     *  \param [in] no input parameters needed
     *
     *   \return [std::vector<short>]  integer ambiguities
     */
    std::vector<short> get_all_integerAmbiguities() const { return _all_integerAmbiguities; };

  /**
     *  \b Description: \n
     * returns a vector containing the multi-band delays (same order as in ncfile)
     * if no  multi-band delays  could be computed (e.g. missing observatioon due to useflag)
     * the value is 666 
     * make sure that fill_all_delay_and_integer_matrix is called befor using this getter 
     * 
     *  \param [in] no input parameters needed
     *
     *   \return [std::vector<short>]  multi-band delays
     */
    ivg::Matrix get_all_group_delays() const { return _all_group_delays; };

#if DEBUG
    static unsigned int count;
    unsigned int ID;
#endif


private:   
    /**
     *  \b Description: \n
     * fills the _triangles vector with all found triangles
     * 
     * THE CLASS TRIANGLE INCLUDES POINTERS TO PRIVATE ATTRIBUTES OF THIS !
     * IF A COPY OF THIS IS CREATED and you want to use the Triangles
     * _createTriangles has to run again !  otherwise SEGFAULTS
     * 
     *  \param [in] [bool] only triangles including the reference station are added
     *
     */
    void _createTriangles(bool includesRefStation = true);

    /**
     *  \b Description: \n
     * call this function after baseline.shift has been changed
     * The following functions affect shift
     *    - compute_baselinewise_integer_ambiguities
     *    - apply_closure_condition
     * 
     * updates _baselines.delays, baseline.residulas and baseline.integerAmbiguities
     * in each baseline:
     *      integerAmbiguities += shift 
     *      residulas += shift * spacing
     *      delays += shift * spacing
     * 
     * Afters shifting those variables shift is reseted to zero
     */
    void _shift_residuals_and_delays();

    std::vector<Triangle> _triangles;
    std::map<std::pair<std::string, std::string>, Baseline> _baselines;

    unsigned short _n_bl;
    unsigned short _n_tr;
    std::string _refSta;

    // -those attributes are initalised with fill_all_delay_and_integer_matrix
    //  -> no need to copy the baseline wise saved delays serveral times into one
    //     vector
    // -the full matrix is needed for writing nc files and computing ionosphere
    // -there shall be no constructor doing this
    std::vector<short> _all_integerAmbiguities;
    ivg::Matrix _all_group_delays;
    ivg::Matrix _all_group_delays_sigmas;
    ivg::Matrix _all_group_delays_residuals;
    
    // Attributes copied from session
    unsigned int _nobs_orig;
    ivg::Matrix _effFreq;
    std::string _name;

};

#endif /* NETWORK_H */

