#ifndef SOURCE_ANALYSIS_H
#define	SOURCE_ANALYSIS_H

#include <QWidget>
#include "plot.h"
#include <vector>
#include <iterator>
#include <sstream>
#include <cstdlib>
#include <map>
#include <utility>
#include "tictoc.h"
#include "matrix.h"
#include "trf.h"
#include "fit.h"
#include <libconfig.h++>
#include "param_list.h"
#include "crf.h"

#include "statistics.h"

/**
* @brief Station_analysis class
* @author AI - bakkari developer team
* @date 2015-07-10
* @version 0.1
*/

//using namespace ivg;
using namespace libconfig;
using namespace std;

namespace ivgat
{
   
class Source_analysis {
    
public:
    
    Source_analysis(vector<ivg::Crf*> crfs, vector<ivg::Param_list*> paramlists);
    
    ivg::Matrix get_erp_data(int dataset_nr, ivg::paramtype type, int order);
    
    virtual ~Source_analysis();
    
private:

    vector<ivg::Crf*> _crf_ptrs;
    vector<ivg::Param_list*> _param_list_ptrs;
    
    
};

} // # namespace ivgat

#endif	/* SOURCE_ANALYSIS_H */

