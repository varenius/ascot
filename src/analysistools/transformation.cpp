#include "transformation.h"

namespace ivgat
{

// ...........................................................................
Transformation::Transformation(t_type type) : _type(type)
// ...........................................................................
{
    // check desired type of transformation
    if(_type == t_type::cart_3D || _type == t_type::sph_3D)
    {
        // define default rotation centroid for 3D
        _rot_centroid = ivg::Matrix(3,1,0.0);
        // define default approximation values for the transformation parameter estimation
        _approx = ivg::Matrix(7,1,0.0);
        // in general scale = 1.0
        _approx(6) = 1.0;
    }
    else
        throw runtime_error( "Transformation::Transformation(t_type type): Only cartesian and spheric 3D transformation implemented yet!" );
    
}
// ...........................................................................
void Transformation::set_systems(ivg::Matrix &system1, ivg::Matrix &system2) 
// ...........................................................................
{
    
    if(_type == t_type::cart_3D || _type == t_type::sph_3D)
    {
        // check if dimensions of both systems are correct: same length and 3 columns
        if(system1.rows() == system2.rows() && system1.cols() == system2.cols() && system1.cols() == 3)
        {
            // in case of cartesian coordinates we don't need any prior transformations
            if(_type == t_type::cart_3D)
            {
                _system1 = system1;
                _system2 = system2;
            }
            else if(_type == t_type::sph_3D)
            {
                // in case of spheric coordinates we need to tranform them to cartesian in a first step
                _system1 = transform_sph_to_cart(system1);
                _system2 = transform_sph_to_cart(system2);
            }

            // setting the number of identical points
            _n_corr_pnt = _system1.length();
           _data_useable = true;
        }
        else
            throw runtime_error( "void Transformation::set_systems(ivg::Matrix &system1, ivg::Matrix &system2) : Wrong matrix dimensions of system1 and system2" );
        
    }
    else if(_type == t_type::cart_2D || _type == t_type::pol_2D)
    {
        // not implemented yet
    }
    
}
// ...........................................................................
string Transformation::estimate_parameter(vector<t_param> tp, ivg::Matrix &param, ivg::Matrix &accu, ivg::Matrix &corr) 
// ...........................................................................
{
    if(_data_useable)
    {
        // threshold for aborting the iteration
        #define THRESHOLD 1e-6
        #define MAX_ITER 100

        if(_type == t_type::cart_2D || _type == t_type::pol_2D)
            if(find(tp.begin(),tp.end(),t_param::rz) != tp.end() || find(tp.begin(),tp.end(),t_param::tz) != tp.end())
               throw runtime_error( "ivg::Matrix Transformation::estimate_parameter(vector<t_param> tp, double threshold) : Estimating of rz and/or tz not possible in 2D" );

        // need to check if order is correct (tx,ty,tz,rx,ry,rz,s)
        for(int i=0; i<tp.size()-1; i++)
            if(tp.at(i) > tp.at(i+1))
                throw runtime_error( "ivg::Matrix Transformation::estimate_parameter(vector<t_param> tp, double threshold): Transformation parameters not in ascending order!" );

        // get initialized default matrices
        // roation centroid
        ivg::Matrix rc = _rot_centroid;
        // approximation matrix
        ivg::Matrix approx = _approx;
        // weight matrix
        ivg::Matrix P_ll;
        P_ll.eye(3*_n_corr_pnt);

        int cnt=0;
        while(1)
        {
            // approximation values    
            double x0 = approx(0);
            double y0 = approx(1);
            double z0 = approx(2);
            double ex = approx(3);
            double ey = approx(4);
            double ez = approx(5);
            double m = approx(6);

            // generating design matrix and observation vector
            ivg::Matrix A(3*_n_corr_pnt, 7, 0.0);
            ivg::Matrix l(3*_n_corr_pnt, 1, 0.0);

            for(int i=1; i<=_n_corr_pnt; i++)
            {
                A(i*3-3,0)=-1;
                A(i*3-2,1)=-1;
                A(i*3-1,2)=-1;
                A(i*3-3,3)=-m*((cos(ex)*sin(ey)*cos(ez)-sin(ex)*sin(ez))*(_system1(i-1,1)-rc(1))+(sin(ex)*sin(ey)*cos(ez)+cos(ex)*sin(ey))*(_system1(i-1,2)-rc(2)));
                A(i*3-3,4)=-m*((-sin(ey)*cos(ez))*(_system1(i-1,0)-rc(0))+(sin(ex)*cos(ey)*cos(ez))*(_system1(i-1,1)-rc(1))+(-cos(ex)*cos(ey)*cos(ez))*(_system1(i-1,2)-rc(2)));
                A(i*3-3,5)=-m*((-cos(ey)*sin(ez))*(_system1(i-1,0)-rc(0))+(-sin(ex)*sin(ey)*sin(ez)+cos(ex)*cos(ez))*(_system1(i-1,1)-rc(1))+(+cos(ex)*sin(ey)*sin(ez)+sin(ex)*cos(ex))*(_system1(i-1,2)-rc(2)));
                A(i*3-3,6)=-((cos(ey)*cos(ez))*(_system1(i-1,0)-rc(0))+(sin(ex)*sin(ey)*cos(ez)+cos(ex)*sin(ez))*(_system1(i-1,1)-rc(1))+(-cos(ex)*sin(ey)*cos(ez)+sin(ex)*sin(ez))*(_system1(i-1,2)-rc(2)));
                A(i*3-2,3)=-m*((-cos(ex)*sin(ey)*sin(ez)-sin(ex)*cos(ez))*(_system1(i-1,1)-rc(1))+(-sin(ex)*sin(ey)*sin(ez)+cos(ex)*cos(ez))*(_system1(i-1,2)-rc(2)));
                A(i*3-2,4)=-m*((sin(ey)*sin(ez))*(_system1(i-1,0)-rc(0))+(-sin(ex)*cos(ey)*sin(ez))*(_system1(i-1,1)-rc(1))+(cos(ex)*cos(ey)*sin(ez))*(_system1(i-1,2)-rc(2)));
                A(i*3-2,5)=-m*((-cos(ey)*cos(ez))*(_system1(i-1,0)-rc(0))+(-sin(ex)*sin(ey)*cos(ez)-cos(ex)*sin(ez))*(_system1(i-1,1)-rc(1))+(cos(ex)*sin(ey)*cos(ez)+sin(ex)*sin(ez))*(_system1(i-1,2)-rc(2)));
                A(i*3-2,6)=-((-cos(ey)*sin(ez))*(_system1(i-1,0)-rc(0))+(-sin(ex)*sin(ey)*sin(ez)+cos(ex)*cos(ez))*(_system1(i-1,1)-rc(1))+(cos(ex)*sin(ey)*sin(ez)+sin(ex)*cos(ez))*(_system1(i-1,2)-rc(2)));
                A(i*3-1,3)=-m*((-cos(ex)*cos(ey))*(_system1(i-1,1)-rc(1))+(-sin(ex)*cos(ey))*(_system1(i-1,2)-rc(2)));
                A(i*3-1,4)=-m*((cos(ey))*(_system1(i-1,0)-rc(0))+(-sin(ex)*(-sin(ey)))*(_system1(i-1,1)-rc(1))+(-cos(ex)*sin(ey))*(_system1(i-1,2)-rc(2)));
                A(i*3-1,5)=0;
                A(i*3-1,6)=-((sin(ey))*(_system1(i-1,0)-rc(0))+(-sin(ex)*cos(ey))*(_system1(i-1,1)-rc(1))+(cos(ex)*cos(ey))*(_system1(i-1,2)-rc(2)));

                l(i*3-3,0)=-rc(0)+_system2(i-1,0)-x0-m*((cos(ey)*cos(ez))*(_system1(i-1,0)-rc(0))+(sin(ex)*sin(ey)*cos(ez)+cos(ex)*sin(ez))*(_system1(i-1,1)-rc(1))+(-cos(ex)*sin(ey)*cos(ez)+sin(ex)*sin(ez))*(_system1(i-1,2)-rc(2)));
                l(i*3-2,0)=-rc(1)+_system2(i-1,1)-y0-m*((-cos(ey)*sin(ez))*(_system1(i-1,0)-rc(0))+(-sin(ex)*sin(ey)*sin(ez)+cos(ex)*cos(ez))*(_system1(i-1,1)-rc(1))+(cos(ex)*sin(ey)*sin(ez)+sin(ex)*cos(ez))*(_system1(i-1,2)-rc(2)));
                l(i*3-1,0)=-rc(2)+_system2(i-1,2)-z0-m*((sin(ey))*(_system1(i-1,0)-rc(0))+(-sin(ex)*cos(ey))*(_system1(i-1,1)-rc(1))+(cos(ex)*cos(ey))*(_system1(i-1,2)-rc(2)));
            }

            // only take desired paramaters
            vector<int> indexes;
            for(int p=0; p<t_param::MAXTPARAM; p++)
                if(find(tp.begin(),tp.end(),p) == tp.end())
                    indexes.push_back(p);

            A.rem_c(indexes);

            // least-square-adjustment
            l = l*(-1.0);
            ivg::Lsa lsa(A,l,P_ll);
            lsa.solve_neq(0);
            ivg::Matrix delta_x = lsa.get_parameters();

            // only update estimated parameters
            for(int p=0; p<tp.size(); p++)
                approx((int)tp.at(p)) += delta_x(p);

            // determinte check-variable
            ivg::Matrix test = delta_x.abs().sum_col() / delta_x.length();

            cnt++;
            // if addition small enough, return results
            if(test(0) < THRESHOLD)
            {
                ivg::Matrix results(tp.size(),1,0.0);
                for(int p=0; p<tp.size(); p++)
                    results(p) = approx((int)tp.at(p));

                param = results;
                accu = lsa.get_vcm().diag().sqrt();
                corr = lsa.get_correlation_matrix();
                
                return get_info_block(tp,param,accu,corr);
            }
            else if(cnt > MAX_ITER)
                throw runtime_error( "ivg::Matrix Transformation::transform(vector<t_param> tp, double threshold): Estimation of transformation parameters not converging");

        }
    }
    
    return "";
}
// ...........................................................................
ivg::Matrix Transformation::get_residuals(vector<t_param> tp_enum, ivg::Matrix tp_passed)
// ...........................................................................
{
    // need to check if order is correct (tx,ty,tz,rx,ry,rz,s)
    for(int i=0; i<tp_enum.size()-1; i++)
        if(tp_enum.at(i) > tp_enum.at(i+1))
            throw runtime_error( "ivg::Matrix Transformation::get_residuals(vector<t_param> tp_enum, ivg::Matrix tp_passed): Transformation parameters not in ascending order!" );
    
    // roation centroid
    ivg::Matrix rc = _rot_centroid;
    
    // transformation parameters not passed need to be set to initial _approx values
    ivg::Matrix tp = _approx;
    for(int p=0; p<tp_enum.size(); p++)
        tp((int)tp_enum.at(p)) = tp_passed(p);
        
    // rotated system 1
    ivg::Matrix rot_system1(_n_corr_pnt,3,0.0);
        
    for(int i=1; i<=_n_corr_pnt; i++)
    {
        rot_system1(i-1,1)=rc(2)+tp(1)+tp(6)*((-cos(tp(4))*sin(tp(5)))*(_system1(i-1,0)-rc(0))+(-sin(tp(3))*sin(tp(4))*sin(tp(5))+cos(tp(3))*cos(tp(5)))*(_system1(i-1,1)-rc(1))+(cos(tp(3))*sin(tp(4))*sin(tp(5))+sin(tp(3))*cos(tp(5)))*(_system1(i-1,2)-rc(2)));
        rot_system1(i-1,0)=rc(0)+tp(0)+tp(6)*((cos(tp(4))*cos(tp(5)))*(_system1(i-1,0)-rc(0))+(sin(tp(3))*sin(tp(4))*cos(tp(5))+cos(tp(3))*sin(tp(5)))*(_system1(i-1,1)-rc(1))+(-cos(tp(3))*sin(tp(4))*cos(tp(5))+sin(tp(3))*sin(tp(5)))*(_system1(i-1,2)-rc(2)));
        rot_system1(i-1,2)=rc(2)+tp(2)+tp(6)*((sin(tp(4)))*(_system1(i-1,0)-rc(0))+(-sin(tp(3))*cos(tp(4)))*(_system1(i-1,1)-rc(1))+(cos(tp(3))*cos(tp(4)))*(_system1(i-1,2)-rc(2)));
    }
    
    ivg::Matrix resid = _system2 - rot_system1;
    return resid;
    
}
// ...........................................................................
ivg::Matrix Transformation::transform_sph_to_cart(ivg::Matrix &sph)
// ...........................................................................
{
    // check if dimensions match the expectation [nx3]
    if(sph.cols() == 3)
    {
        ivg::Matrix xyz(sph.rows(),3,0.0);

        for(int i=0; i<sph.rows(); i++)
        {
            double azi = sph(i,0);
            double ele = sph(i,1);
            double r = sph(i,2);
            
            xyz(i,0) = r * cos(ele) * cos(azi);
            xyz(i,1) = r * cos(ele) * sin(azi);
            xyz(i,2) = r * sin(ele);
        }
    
        return xyz;
    }
    else
        throw runtime_error( "ivg::Matrix Transformation::transform_sph_to_cart(ivg::Matrix &sph): Wrong dimension of matrix sph. [nx3] expected.");
}
// ...........................................................................
ivg::Matrix Transformation::transform_cart_to_sph(ivg::Matrix &xyz)
// ...........................................................................
{
    // check if dimensions match the expectation [nx3]
    if(xyz.cols() == 3)
    {
        ivg::Matrix sph(xyz.rows(),3,0.0);

        for(int i=0; i<xyz.rows(); i++)
        {
            double x = xyz(i,0);
            double y = xyz(i,1);
            double z = xyz(i,2);
            
            double dist_xy = sqrt(pow(abs(x),2)+pow(abs(y),2)); 
            
            sph(i,0) = atan2(y,x); // azimut
            sph(i,1) = atan2(z,dist_xy); // elevation
            sph(i,2) = sqrt(pow(abs(dist_xy),2)+pow(abs(z),2)); // radius
        }
    
        return sph;
    }
    else
        throw runtime_error( "ivg::Matrix Transformation::transform_cart_to_sph(ivg::Matrix &xyz): Wrong dimension of matrix xyz. [nx3] expected.");
}
// ...........................................................................
void Transformation::rotate_3d(ivg::Matrix &xyz, double alpha, double beta, double gamma) 
// ...........................................................................
{
    ivg::Matrix rot_xyz,rot_x,rot_y,rot_z;
    rot_x.rot3D_x(alpha);
    rot_y.rot3D_y(beta);
    rot_z.rot3D_z(gamma);
    
    rot_xyz = rot_x*rot_y*rot_z;
        
    for(int i=0; i<xyz.rows(); i++)
        xyz.set_sub(i,0,(rot_xyz*(xyz.get_sub(i,0,i,2).transpose())).transpose());
    
}
// ...........................................................................
string Transformation::get_info_block(vector<t_param> trans_param, ivg::Matrix &param, ivg::Matrix &accu, ivg::Matrix &corr)
// ...........................................................................
{
        stringstream trans_info;
        trans_info << "------------------------------------------------------------------" << endl;
        trans_info << "Estimated Transformation Parameter" << endl;
        trans_info << "------------------------------------------------------------------" << endl;
        double fak=1;
        string unit = "-";
        int cnt=0;
        for(auto &tp: trans_param)
        {
            if( tp <= 2 )
            {
                fak = 1000.0;
                unit = "mm";
            }
            else if( tp > 2 && tp <= 5)
            {
                fak = ivg::rad2mas*1000;
                unit = "microas";
            }
            else if( tp == 6 )
            {
                fak = 1.0;
                unit = "-";
            }
            
            trans_info << ivgat::t_param_name.at(tp) << " ";
            trans_info << right << setw(18) << fixed << setprecision(12) << param(cnt)*fak;
            trans_info << "  +-  ";
            trans_info << setw(10) << left << setprecision(5) << fixed << accu(cnt)*fak << " [" << unit << "]" << endl;
            
            cnt++;
        }
        
        trans_info << "------------------------------------------------------------------" << endl;
        trans_info << "Correlation Matrix" << endl;
        trans_info << "------------------------------------------------------------------" << endl;
        for(int i=0; i<corr.rows(); i++)
        {
            for(int j=0; j<corr.cols(); j++)
                trans_info << setw(5) << right << setprecision(2) << fixed << corr(i,j) << " ";

            trans_info << endl;
        }
        trans_info << "------------------------------------------------------------------" << endl;
    
        
    return trans_info.str();
}

} // namespace ivgat end
