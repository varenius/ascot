#ifndef PARAM_H
#define	PARAM_H

#include <vector>
#include "logger.h"
#include "ivg_const.h"
#include "date.h"
#include "auxfunc.h"
#include "structs.h"
#include <iterator>
#include <iostream>
#include <algorithm>
#include <fstream>
#include <sstream>
#include "tictoc.h"

/**
*
* @brief class Param - base class
* @author AI - bakkari developer team
* @date 2015-03-24
* @version 0.1
*/

using namespace std;
namespace ivg
{

//Parameter Type
                                enum paramtype {        stax,        stay,        staz,    clo,         zwd,         ngr,         egr,          xpo,          ypo,            ut1,         nutx,         nuty,      ra,     dec,     cbr,    nutln,     nutob,    blcl, atbr, MAXPARAM };
const std::vector<std::string> paramtype_str = {      "stax",      "stay",      "staz",  "clo",    "zwd"   ,       "ngr",       "egr",        "xpo",        "ypo",          "ut1",       "nutx",       "nuty",    "ra",   "dec",   "cbr",  "nutln",   "nutob",  "blcl",  "atbr"};
const std::vector<std::string> paramtype_snx(  {      "STAX",      "STAY",      "STAZ",  "CLO",    "TROTOT",    "TGNTOT",    "TGETOT",        "XPO",        "YPO",          "UT1",      "NUT_X",      "NUT_Y", "RS_RA", "RS_DE", "CL_BR", "NUT_LN",  "NUT_OB", "BL_CL", "AT_CR"} );
const std::vector<std::string> paramtype_unit( {      "m"   ,         "m",         "m",    "s",         "m",         "m",         "m",        "mas",        "mas",           "ms",        "mas",        "mas",   "rad",   "rad",     "s",    "mas",     "mas",     "s",     "m"} );
const std::vector<double> param_unit_fac(      { ivg::c*1e-2, ivg::c*1e-2, ivg::c*1e-2, 1.0e-1, ivg::c*1e-1, ivg::c*1e-2, ivg::c*1e-2, ivg::rad2mas, ivg::rad2mas, ivg::rad2s*1e3, ivg::rad2mas, ivg::rad2mas,     1.0,     1.0,  1.0e-1,      1.0,       1.0,  1.0e-1, ivg::c*1e-1} );

//enum paramtype { stax, stay, staz, clo, zwd, ngr, egr, xpo, xpor, ypo, ypor, ut1, lod, nutx, nuty, ra, dec };
//const std::vector<std::string> paramtype_str = { "stax", "stay", "staz", "clo", "zwd", "ngr", "egr", "xpo", "xpor", "ypo", "ypor", "ut1", "lod", "nutx", "nuty", "ra", "dec" };
//const std::vector<std::string> paramtype_snx( { "STAX", "STAY", "STAZ", "CLO", "TROTOT", "TGNTOT", "TGETOT", "XPO", "XPOR", "YPO", "YPOR", "UT1", "LOD", "NUT_X", "NUT_Y", "RS_RA", "RS_DE" } );
//const std::vector<std::string> paramtype_unit( { "m", "m", "m", "-", "m", "m", "m", "mas", "masD", "mas", "masD", "ms", "ms", "mas", "mas", "rad", "rad"} );

class Param
{

    public:

        // ==============================================
        // =============== Constructors: ================
        // ==============================================
        /**
        *  \b Description: \n
        *        Default constructor
        *  \param [in] no input parameters needed
        *  \return An instance of the class 'Param'
        */
        Param();

        /**
        *  \b Description: \n
        *        constructor using two input parameter: parameter type and name
        *  \param [in] [ivg::paramtype] parameter type
        *              [std::string] name (e.g., station or source name)
        *  \return An instance of the class 'Param'
        */
        Param(ivg::paramtype type, string name);

//    Param(ivg::paramtype type, string name, ivg::Date epoch, int order = 0);

        /**
        *  \b Description: \n
        *        constructor using five input parameter:
        *        parameter type and name, epoch, a priori value and degree of polynom
        *  \param [in] [ivg::paramtype] parameter type
        *              [std::string] name (e.g., station or source name)
        *              [ivg::Date] epoch
        *              [double] a priori value
        *              [int] degree of polynom
        *  \return An instance of the class 'Param'
        */
        Param(ivg::paramtype type, string name, ivg::Date epoch, double apriori = 0.0, int order = 0);

        /**
        *  \b Description: \n
        *        constructor using two input parameter: line of an external file (e.g., 
        *        the intern ivg::ASCOT result file) and a boolian to specify whether 
        *        the file contains a priori information or not.
        *         
        *  \param [in] [ivg::paramtype] parameter type
        *              [std::string] line of an external file
        *              [bool] specify whether the file contains a priori information or not 
        *  \return An instance of the class 'Param'
        */        
        Param(const string line, bool apriori);
        
        /**
        *  \b Description: \n
        *        copy constructor
        *  \param [in] [ivg::Param] other 'Param' object
        *  \return An instance of the class 'Param'
        */
        Param(const Param& orig);

        /**
        *  \b Description: \n
        *        Default deconstructor
        *  \param [in] no input parameters needed
        */
        virtual ~Param();


        // ==============================================
        // =========== MEMBER-functions: ================
        // ==============================================
        /**
        *  \b Description: \n
        *        Method to check wheather a test Param is equal to this Param
        *        Please note: a parameter is equal to another parameter, if
        *        the name, the type and the polynomial degree are equal
        *  \param [in] [ivg::Session] test session
        *  \return [bool] true: same session
        *                 false: other session
        */
        bool operator==( const Param test ) const;
        
        void compare_to( Param &test, bool &name, bool &type, bool &order, bool &epoch );

        /**
        *  \b Description: \n
        *        Method to set the epoch
        *  \param [in] [ivg::Date] new epoch
        */
        void set_epoch( ivg::Date epoch )
        {
            _epoch = epoch;
        };
        
        /**
        *  \b Description: \n
        *        Method to set the name of the parameter (e.g. TIGOCONC)
        *  \param [in] [string] new name
        */
        void set_name( string name )
        {
            _name = name;
        };

        /**
        *  \b Description: \n
        *        Method to set the polynomial degree
        *  \param [in] [int] degree of the polynom
        */
        void set_order( int o )
        {
            _order=o;
        };

        /**
        *  \b Description: \n
        *        Method to set an a priori value
        *  \param [in] [double] a priori value
        */
        void set_apriori( double a )
        {
	  
	  _apriori=a;
        };

        /**
        *  \b Description: \n
        *        Method to set a estimate value
        *  \param [in] [double] estimate value
        */
        void set_estimate( double e )
        {
            _estimate = e;
        };

        /**
        *  \b Description: \n
        *        Method to set a estimate value
        *  \param [in] [double] estimate value
        */
        void set_standard_deviation( double std )
        {
            _std = std;
        };

        /**
        *  \b Description: \n
        *        Method to set the standard deviation of the offset constraints
        *  \param [in] [double] standard deviation of the offset constraints
        */
        void set_offset_cnstr_sigma( double c )
        {
            _offset_cnstr_sigma = c;
        };

        /**
        *  \b Description: \n
        *        Method to set the standard deviation of the rate constraints
        *  \param [in] [double] standard deviation of the rate constraints
        */
        void set_rate_cnstr_sigma( double c )
        {
            _rate_cnstr_sigma = c;
        };

        /**
        *  \b Description: \n
        *        Method to change the parameter type
        *  \param [in] [ivg::paramtype] new parameter type
        *  \return [ivg::paramtype] parameter type
        */
        void set_type( ivg::paramtype new_type )
        {
            _type = new_type;
        };
        
        void increment_stacked()
        {
            _stacked += 1;
        }
        
        void set_stacked( int stacked )
        {
            _stacked = stacked;
        }

        /**
        *  \b Description: \n
        *        Method to get the parameter type
        *  \param [in] no input parameters needed
        *  \return [ivg::paramtype] parameter type
        */
        ivg::paramtype get_type()
        {
            return _type;
        };

        /**
        *  \b Description: \n
        *        Method to get the polynomial degree
        *  \param [in] no input parameters needed
        *  \return [int] polynomial degree
        */
        int get_order()
        {
            return _order;
        };
        
        /**
        *  \b Description: \n
        *        Method to get the number how often parameter has been stacked
        *  \param [in] no input parameters needed
        *  \return [int] amount of stacking
        */
        int get_stacked()
        {
            return _stacked;
        };

        /**
        *  \b Description: \n
        *        Method to get the parameter name (e.g., station name or source name)
        *  \param [in] no input parameters needed
        *  \return [std::string] parameter name
        */
        std::string get_name()
        {
            return _name;
        };
        
        /**
        *  \b Description: \n
        *        Method to get the parameter typename (e.g., stax, ut1, dec)
        *  \param [in] no input parameters needed
        *  \return [std::string] typename
        */
        std::string get_typename()
        {
            return paramtype_str.at(_type);
        };
        
        /**
        *  \b Description: \n
        *        Method to get the common parameter unit
        *  \param [in] no input parameters needed
        *  \return [std::string] unit
        */
        std::string get_unit()
        {
            return paramtype_unit.at(_type);
        };

        /**
        *  \b Description: \n
        *        Method to get the a priori value of the parameter
        *  \param [in] no input parameters needed
        *  \return [double] a priori value of the parameter
        */
        double get_apriori()
        {
            return _apriori;
        };

        /**
        *  \b Description: \n
        *        Method to get the estimate
        *  \param [in] no input parameters needed
        *  \return [double] estimate
        */
        double get_estimate()
        {
            return _estimate;
        };

        /**
        *  \b Description: \n
        *        Method to get the standard deviation of the estimate
        *  \param [in] no input parameters needed
        *  \return [double] standard deviation of the estimate
        */
        double get_standard_deviation()
        {
            return _std;
        };

        /**
        *  \b Description: \n
        *        Method to get the epoch of the parameter
        *  \param [in] no input parameters needed
        *  \return [ivg::Date] epoch of parameter
        */
        ivg::Date get_epoch()
        {
            return _epoch;
        };

        /**
        *  \b Description: \n
        *        Method to get the standard deviation of the rate constraints
        *  \param [in] no input parameters needed
        *  \return [double] standard deviation of the rate constraints
        */
        double get_rate_cnstr_sigma()
        {
            return _rate_cnstr_sigma;
        };

        /**
        *  \b Description: \n
        *        Method to get the standard deviation of the offset constraints
        *  \param [in] no input parameters needed
        *  \return [double] standard deviation of the offset constraints
        */
        double get_offset_cnstr_sigma()
        {
            return _offset_cnstr_sigma;
        };

        /**
        *  \b Description: \n
        *        Method to get the standard deviation of the offset and rate constraints
        *  \param [in] no input parameters needed
        *  \param [out] [double] standard deviation of the offset constraints
        *               [double] standard deviation of the rate constraints
        */
        void get_cnstr_sigmas( double & offset, double & rate );


        // to-do: still in use?

//   void set_rate_cnstr_ptr( Param * param_ptr, std::string pos );
//   Param * next_param(){ return _after_para; };
//   Param * prev_param(){ return _previ_para; };

        /**
        *  \b Description: \n
        *        Method to show the parameter
        *  \param [in] no input parameters needed
        */
        void show();

        /**
        *  \b Description: \n
        *        Method to get one line including results of the estimation process for a parameter.
        *        This includes:
        *        (a) index
        *        (b) name
        *        (c) type
        *        (d) polynomial degree
        *        (e) epoch (MJD)
        *        (f) a priori value
        *        (g) estimate
        *        (h) standard deviation
        *  \param [in] no input parameters needed
        */
        string get_resultline(bool apriori);
        void parse_resultline(const string rl, bool apriori);

        /**
        *  \b Description: \n
        *        Method to detect if parameter is from a definied type and 
        *        a defined order
        *  \param [in] types, which the parameter should be checked for
        *  \param [in] orders, which the parameter should be checked for
        *  \return [out] is true, if parameter is one of the defined types/orders
        */
        bool is_type(vector<ivg::paramtype> types, vector<int> orders);

        bool is_type_name(vector<ivg::paramtype> types, vector<std::string> names);
        
        void set_reduce_flag( bool f ){ _reduce = f;};
        bool get_reduce_flag(){ return _reduce; };

    private:

        // ==============================================
        // ====== class variables / attributes:==========
        // ==============================================

        // parameter type, name and polynomial degree
        ivg::paramtype _type;
        string _name;
        int _order;
        int _stacked;

        // a priori value, estimation value, standard deviation (polynom),
        // standard deviation of the offset constraint and epoch of a parameter
        double _apriori;
        double _estimate;
        double _std;
        double _offset_cnstr_sigma;
        ivg::Date _epoch;

        // CPWLF/B-spline parameters
        Param * _previ_para;
        Param * _after_para;
        double _rate_cnstr_sigma;

        bool _reduce;

};

}

#endif	/* PARAM_H */

