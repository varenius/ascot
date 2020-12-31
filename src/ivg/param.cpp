#include "param.h"

namespace ivg
{
// ...........................................................................
Param::Param()
{
}
// ...........................................................................
Param::Param(ivg::paramtype type, string name)
// ...........................................................................
{
    _name = name;
    _type = type;
    _offset_cnstr_sigma = 0.0;
    _rate_cnstr_sigma = 0.0;
    _estimate = 0.0;
    _std = 0.0;
    _order = 0;
    _stacked = 1;
    _previ_para = NULL;
    _after_para = NULL;
    _reduce = false;
}
// ...........................................................................
Param::Param(ivg::paramtype type, string name, ivg::Date epoch, double apriori, 
             int order) : Param(type,name)
// ...........................................................................
{
    _epoch = epoch;
    _apriori = apriori;
    _order = order;
}
// ...........................................................................
Param::Param(const string line, bool apriori)
// ...........................................................................
{
    std::istringstream in(line);

    string name, paramtype;
    int order;
    double epoch, apri_plus_esti_times_factor, esti_times_factor, apri_times_factor, std_times_factor;
    
    if (apriori) 
    {
        in >> _stacked >> name >> paramtype >> order >> epoch >> apri_plus_esti_times_factor >> apri_times_factor >> esti_times_factor >> std_times_factor;
    } 
    else 
    {
        in >> name >> paramtype >> order >> epoch >> apri_plus_esti_times_factor >> esti_times_factor >> std_times_factor;
        apri_times_factor = apri_plus_esti_times_factor - esti_times_factor;
    }

    vector<string>::const_iterator pos_iter = find(paramtype_str.begin(), paramtype_str.end(), paramtype);
    int pos = pos_iter - paramtype_str.begin();

    *this = Param((ivg::paramtype)pos, name, ivg::Date(epoch), apri_times_factor, order);
    _estimate = esti_times_factor;
    _std = std_times_factor;
}
// ...........................................................................
Param::Param(const Param& orig)
// ...........................................................................
{
    _name = orig._name;
    _type = orig._type;
    _epoch = orig._epoch;
    _apriori = orig._apriori;
    _estimate = orig._estimate;
    _std = orig._std;
    _order = orig._order;
    _stacked = orig._stacked;
    _offset_cnstr_sigma = orig._offset_cnstr_sigma;
    _rate_cnstr_sigma = orig._rate_cnstr_sigma;
    _previ_para = orig._previ_para;
    _after_para = orig._after_para;
    _reduce = orig._reduce;
}
// ...........................................................................
Param::~Param()
// ...........................................................................
{}
// ...........................................................................
bool Param::operator==( const Param test ) const
// ...........................................................................
{
    bool out = false;
    if( _name == test._name && _type == test._type && _order == test._order )
        out = true;

    return out;
}
// ...........................................................................
void Param::compare_to( Param &test, bool &name, bool &type, bool &order, 
                        bool &epoch )
// ...........................................................................
{
    name = false;
    type = false;
    order = false;
    epoch = false;
    
    if( _name == test._name)
        name = true;
    if(_type == test._type)
        type = true;
    if(_order == test._order )
        order = true;
    if(_epoch == test._epoch)
        epoch = true;
    
}
// ...........................................................................
void Param::show()
// ...........................................................................
{
    cout << get_resultline(true);
    cout << endl;
}
// ...........................................................................
string Param::get_resultline(bool apriori)
// ...........................................................................
{

    ostringstream ss;
    if(apriori)
    ss << setw(4) << _stacked << " ";
    ss << setfill(' ') << setw(8) << left << _name << " "
        << setfill(' ') << setw(5) << left << paramtype_str.at(_type) << " "
        << setfill(' ') << setw(2) << left << _order << " "
        << setfill(' ') << setw(10) << left << setprecision(5) << fixed 
        << _epoch.get_double_mjd() << " "
        << setfill(' ') << setw(22) << right << setprecision(14) << scientific 
        << (_apriori+_estimate)*ivg::param_unit_fac.at(_type) << " ";
    if(apriori)
        ss << setfill(' ') << setw(22) << right << setprecision(14) <<scientific 
           << _apriori*ivg::param_unit_fac.at(_type) << " ";
    ss  << setfill(' ') << setw(22) << right << setprecision(14) << scientific 
        << _estimate*ivg::param_unit_fac.at(_type) << " "
        << setfill(' ') << setw(17) << right << setprecision(11) << scientific 
        << _std*ivg::param_unit_fac.at(_type) << " "
        << setw(3) << left << ivg::paramtype_unit.at(_type) << " ";
//        << setfill(' ') << setw(22) << right << setprecision(14) << scientific << _apriori << " ";
      
    return ss.str();
}
/*
// ...........................................................................
void Param::set_rate_cnstr_ptr( Param * param_ptr, std::string pos )
// ...........................................................................
{
   cerr << "/// " << param_ptr << " <= " << this << endl;
	if( pos == "previous" )
		_previ_para = param_ptr;
	else if ( pos == "after" )
		_after_para = param_ptr;
	else
	{
		stringstream errormessage;
		errormessage << "Param::set_rate_cnstr_ptr( Param * param_ptr, std::string pos ) "
					 << "wrong input parameter! "
					 << " Exiting";
		throw runtime_error( errormessage.str() );
	}
}
*/
// ...........................................................................
void Param::get_cnstr_sigmas( double & offset, double & rate )
// ...........................................................................
{
    offset = _offset_cnstr_sigma;
    rate   = _rate_cnstr_sigma;
}
// ...........................................................................
bool Param::is_type(vector<ivg::paramtype> types, vector<int> orders)
// ...........................................................................
{    
    for(int i=0; i<types.size(); i++)
    {        
        for(int j=0; j<orders.size(); j++)
        {
            if(_type == types.at(i) && _order == orders.at(j))
                return true;
        }
    }    
    return false;
}

// ...........................................................................
bool Param::is_type_name(vector<ivg::paramtype> types, vector<std::string> names)
// ...........................................................................
{
    
    for(int i=0; i<types.size(); i++)
    {
        
        for(int j=0; j<names.size(); j++)
        {
            if(_type == types.at(i) && _name == names.at(j))
                return true;
        }
    }
    
    return false;
}


}
