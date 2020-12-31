#ifndef STATIONS_H
#define STATIONS_H

#include "matrix.h"
#include "date.h"
#include <map>
#include <iostream>
#include <string>

/**
*
* @brief abstract class Station - base class
* @author SH - bakkari developer team
* @date 2015-03-24
* @version 0.1
*
*/

using namespace std;

namespace ivg
{

//Station-Name-Type
enum staname { MINSTA, domes_no, cdp , ivs_name, lettercode, description, corres, ant_name, equip_id, MAXSTA};

// ===========================================================================
class Station
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
        *  \return An instance of the class 'Station'
        */
        Station();

        /**
        *  s\b Description: \n
        *        Constructor using four parameters:
        *        position and velocity matrix, reference epoch and names
        *  \param [in] [ivg::Matrix] position vector
        *               [ivg::Matrix] velocity vector
        *               [ivg::Matrix] reference epoch
        *               [std::map<std::string, std::string>] station names
        *               (possible types: name, domes no, cdp, ivs name, lettercode, description)
        *  \return An instance of the class 'Station'
        */
        Station( ivg::Matrix xyz0, ivg::Matrix vel0, ivg::Date refepoch,
                 std::vector<ivg::Date> breakepochs,
                 std::map<ivg::staname, std::string> names );

        /**
        *  \b Description: \n
        *        Default deconstructor
        *  \param [in] no input parameters needed
        */
        ~Station();


        // ==============================================
        // =========== MEMBER-functions: ================
        // ==============================================

        // getter
        /**
        *  \b Description: \n
        *        Method to get position matrix
        *  \param [in] no input parameters needed
        *  \return [Matrix] position vector
        */
        ivg::Matrix get_xyz0() const;

        /**
        *  \b Description: \n
        *        Method to get the standard deviation of the initial XYZ station coordinations
        *  \param [in] no input parameters needed
        *  \return [Matrix] position vector
        */
        ivg::Matrix get_xyz0_std() const;        
        
        /**
        *  \b Description: \n
        *        Method to get velocity matrix
        *  \param [in] no input parameters needed
        *  \return [Matrix] velocity matrix
        */
        ivg::Matrix get_vel0() const;

        /**
        *  \b Description: \n
        *        Method to get reference epoch
        *  \param [in] no input parameters needed
        *  \return [Matrix] reference epoch
        */
        ivg::Date get_refepoch() const;

        /**
        *  \b Description: \n
        *        Method to get station name of selected name type
        *  \param [in] no input parameters needed
        *  \return [std::string] station name of selected name type
        */
        string get_name( ivg::staname type = MAXSTA);

        /**
        *  \b Description: \n
        *        Method to get names
        *        (name, domes no, cdp, ivs name, lettercode, description)
        *  \param [in] no input parameters needed
        *  \return [std::map<std::string, std::string>] station names of all name types
        */
        std::map<ivg::staname, std::string> get_names( );

        /**
        *  \b Description: \n
        *        Method to list names
        *        (name, domes no, cdp, ivs name, lettercode, description)
        *  \param [in] no input parameters needed
        */
        void show_names( );

        /**
        *  \b Description: \n
        *        Method to modify position matrix
        *  \param [in] [Matrix] new position matrix
        */
        void set_xyz0( Matrix &xyz );

        /**
        *  \b Description: \n
        *        Method to modify velocity matrix
        *  \param [in] [Matrix] new velocity matrix
        */
        void set_vel0( Matrix &vel );

        /**
        *  \b Description: \n
        *        Method to modify reference epoch
        *  \param [in] [Matrix] new reference epoch
        */
        void set_refepoch( ivg::Date refepoch );

        /**
        *  \b Description: \n
        *        Method to modify station name of selected name type
        *  \param [in] [std::string] new station name of selected name type
        */
        void set_name(ivg::staname type, string value);

        /**
        *  \b Description: \n
        *        Pure virtual method to show station
        *  \param [in] no input parameters needed
        */
        virtual void show() = 0;

        /**
        *  \b Description: \n
        *        Pure virtual method to calculate position vector
        *  \param [in] [Matrix] epoch of interest
        *  \return [Matrix] [3 x 1] position vector
        */
        virtual ivg::Matrix calc_xyz( ivg::Date epoch ) = 0;

        /**
        *  \b Description: \n
        *        Pure virtual method to get elevation angle (el) and azimuth (az)
        *  \param [in] no input parameters needed
        *  \return [Matrix] [n x 2] matrix with el-az
        */
//	virtual ivg::Matrix el_az() = 0;



    protected:

        // ==============================================
        // ======== class variables / attributes: =======
        // ==============================================

        // station coordinate (position) matrix
        ivg::Matrix _pos0;
        ivg::Matrix _pos0_std;
        ivg::Matrix _llh0;
  
        // velocity matrix
        ivg::Matrix _vel0;
        ivg::Matrix _vel0_std;

        // reference epochs
        ivg::Date _refepoch;

        // break epochs
        std::vector<ivg::Date> _discontinuity;

        // station name(s)
        map<ivg::staname,string> _names;


};

} // # namespace ivg

#endif // STATIONS_H
