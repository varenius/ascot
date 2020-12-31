#include "source_analysis.h"

namespace ivgat
{
   
Source_analysis::Source_analysis(vector<ivg::Crf*> crfs, vector<ivg::Param_list*> paramlists)
{
    // backup comments
    //    ivg::Session gsf2014a, opa2015a;
    //    snxizer.load( &gsf2014a, &setup, dir+"/snx/191_12_55_47.snx"); // opa2015a 
    //    snxizer.load( &opa2015a, &setup, dir+"/snx/191_12_57_00.snx"); // gsf2014a

    // get pointer to crf to generate vector of crfs
    //    vector<ivg::Crf*> crfs = {gsf2014a.get_crf_ptr(), opa2015a.get_crf_ptr()};
    //    vector<ivg::Param_list*> paramlists = {gsf2014a.get_param_list_ptr(), opa2015a.get_param_list_ptr()};

    //    ivgat::Source_analysis sourcernizer(crfs,paramlists);
    
    
    if(crfs.size() == paramlists.size())
    {
        _crf_ptrs = crfs;
        _param_list_ptrs = paramlists;
        
    }
    else
        throw runtime_error( "Source_analysis::Source_analysis(vector<ivg::Crf*> crfs, vector<ivg::Param_list*> paramlists): Unequal length of crfs and paramlists." );
    
}

ivg::Matrix Source_analysis::get_erp_data(int dataset_nr, ivg::paramtype type, int order) {
    
    ivg::Param_list *tmp_list = _param_list_ptrs.at(dataset_nr);
    
    ivg::Matrix plot_data;
    for(std::vector<ivg::Param>::iterator param = tmp_list->begin(); param != tmp_list->end(); ++param)
    {
        if(param->is_type({type}, {order}))
        {
            ivg::Matrix row(1,3,0.0);
            row(0,0) = param->get_epoch().get_double_mjd();
            row(0,1) = param->get_estimate();
            row(0,2) = param->get_apriori();
            plot_data.append_rows(row);
        } 
    }
    
    
    return plot_data;
}

Source_analysis::~Source_analysis() {
}

 
}