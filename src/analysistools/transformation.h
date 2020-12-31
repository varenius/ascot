/** Transformation Class  
 *   
 *  @author Andreas Iddink 
 *  @date 22 Sep 2015
 * 
 */

#ifndef TRANSFORMATION_H
#define	TRANSFORMATION_H

#include "ivg_const.h"
#include "matrix.h"
#include "lsa.h"

namespace ivgat
{

/// Implemented Transformation Types
enum t_type {cart_2D, cart_3D, pol_2D, sph_3D, MAXTTYPE};

/// Implemented Transformation Parameter s
enum t_param {tx,ty,tz,rx,ry,rz,s,MAXTPARAM};
const vector<string> t_param_name = {"tx","ty","tz","rx","ry","rz","s"};
const vector<string> t_param_unit = {"m","m","m","rad","rad","rad","-"};

class Transformation {
    
public:
    
    /**
    *  \b Description: \n
    *        Constructor based on Enumerator Transformation Type.
    *
    *  \param [in] type - e.g. ivgat::t_type::cart_3D 
    */
    Transformation(t_type type);
    
    /**
    *  \b Description: \n
    *        Method to set the right-handed 3D cartesian coordinates. Two matrices (n x 3)
    *        with coordinates in the origin (x y z).
    *
    *  \param [in] system1 - nx3 or nx2, depends on t_type, see constructor
    *  \param [in] system2 - nx3 or nx2, depends on t_type, see constructor
    */
    void set_systems(ivg::Matrix &system1, ivg::Matrix &system2);
    
    /**
    *  \b Description: \n
    *        Method to estimate the transformation parameter between the with "set_systems" set systems 1 and 2.
    *        Overdetermined cartesian 3D similarity transformation based on iterative Gauss-Markov model
    *
    *  \param [in] vector of type t_param - transformation parameters to estimate
    *  \param [out] param - call by reference - returns estimated parameter (single column)
    *  \param [out] accu - call by reference - returns standard deviation of the estimated parameters (single column) 
    *  \param [out] corr - call by refrence - returns correlation matrix of the estimated parameters 
    */
    string estimate_parameter(vector<t_param> tp, ivg::Matrix &param, ivg::Matrix &accu, ivg::Matrix &corr);
    
    /**
    *  \b Description: \n
    *        Method to get the residuals after transformation both systems on each other.
    *        resids = system2 - f(system1,param) 
    *
    *  \param [in] vector tp_enum, containing parameter types for the transformation
    *  \param [in] parameter matrix, containing the values, defined in vector tp_enum
    *  \return matrix nx 
    */
    ivg::Matrix get_residuals(vector<t_param> tp_enum, ivg::Matrix tp);
    
    /**
    *  \b Description: \n
    *        Method to transform from spherical to cartesian coordinates.
    *
    *  \param [in] sph - matrix - nx3 - with [azimut, elevation, radius]
    *  \return xyz - matrix - nx3 - with transformed [x, y, z]
    */
    ivg::Matrix transform_sph_to_cart(ivg::Matrix &sph);
    
    /**
    *  \b Description: \n
    *        Method to transform from cartesian to spherical coordinates.
    *
    *  \param [in] xyz - matrix - nx3 - with [x, y, z]
    *  \return sph - matrix - nx3 - with transformed [azimut, elevation, radius]
    */
    ivg::Matrix transform_cart_to_sph(ivg::Matrix &xyz);
    
    /**
    *  \b Description: \n
    *        Method to rotate a matrix - nx3 - with given angles alpha(x), beta(y), gamma(z)
    *        Unit of angles RADIAN. Call by reference.
    * 
    *  \param [in] xyz - call by reference - matrix - nx3 - with [x, y, z]
    *  \param [in] alpha, beta, gamma - rotation angles in radian 
    *  \return sph - matrix - nx3 - with transformed [azimut, elevation, radius]
    */
    void rotate_3d(ivg::Matrix &xyz, double alpha, double beta, double gamma);
    
    string get_info_block(vector<t_param> trans_param, ivg::Matrix &param, ivg::Matrix &accu, ivg::Matrix &corr);
    
private:

    // transformation type, e.g. ivgat::t_type::cart_3D 
    t_type _type;
    
    // datasets with equal size [nx3]
    ivg::Matrix _system1,_system2;
    
    // amount of corresponding points
    int _n_corr_pnt;
    
    // variable to check if the systems can be used for operations
    bool _data_useable;

    // default rotation centroid and approximation values for transformation
    ivg::Matrix _rot_centroid;
    ivg::Matrix _approx;
    
    // we have to do this, otherwise compilation error
    ivg::Lsa nonsense;
    
};


} // namespace ivgat end

#endif	/* TRANSFORMATION_H */

