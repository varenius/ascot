#include "ascotizer.h"

Ascotizer::Ascotizer() {
}

Ascotizer::Ascotizer(threadtype type, int session_idx, ivg::Session &session) : _type(type), _idx(session_idx), _session(session)
{

}

void Ascotizer::run(){
    
    // get setup of current sessions

    Setting *setup = _session.get_setup();
    // get sessions_type (vgosdb, ngs, snx)
   
    
    if(_type == threadtype::LOAD)
    {
  
        string session_type = (const char *)get_list_element((*setup)["datadirs"],(*setup)["session_type"])[1];

        ivg::Masterfile masterfile((*setup)["definitions"]["masterfiles"], ivg::mastertype::both);

        ivg::Session_inout sessionizer(session_type,masterfile);

        // hard coded version 4 for NGS
        int year = masterfile.get_session_info(_session.get_name()).date.get_int_year();
            if(year == 0){
                year = stoi(_session.get_name().substr(0,2));
                if (year<79)
                    year += 2000;
                else
                    year += 1900;
            }
	    string vgos_dir = (const char *)get_list_element((*setup)["datadirs"],(*setup)["session_type"])[2];
            vgos_dir += "/" + std::to_string(year) + "/" +_session.get_name()  + "/";
	    ivg::Wrapper wrapper;
            if ( (bool) (*setup) ["use_wrapper"]){
	      wrapper = ivg::Wrapper(vgos_dir,_session.get_name(),(const char *)(*setup)["vgosdb_editing"]);
                sessionizer.setWrapper_ptr(&wrapper);
            }
        sessionizer.load( &_session, setup, _session.get_name(), "-");

    }
    else if(_type == threadtype::INIT)
    {
        _session.init_vgosdb_ngs_solution();
        
        // modify the parameterization
        _session.modify_parameterization();
                
        // reduce parameter and constrain them
        _session.reduce_and_constrain();

        // perform least squares solution
        // contains create_nnr_nnt_equations
        _session.solve();     
        
    }
    else if(_type == threadtype::RESIDUALS)
    {   
        _session.create_solution_info();
    }
    else if(_type == threadtype::EXPORT)
    {
        string session_type = (const char *)get_list_element((*setup)["datadirs"],(*setup)["session_type"])[1];
        ivg::Masterfile masterfile((*setup)["definitions"]["masterfiles"], ivg::mastertype::both);
        ivg::Session_inout sessionizer(session_type,masterfile);
        
        // write each session (solution) results in SNX
        string snx_dir = (*setup)["export_snx"]["dir"];
        string version = (*setup)["export_snx"]["version"];

        if( version == "" )
            sessionizer.write_snx(&_session, snx_dir+"/"+_session.get_name()+".snx" );
        else
            sessionizer.write_snx(&_session, snx_dir+"/"+_session.get_name()+"_"+version+".snx" );
    }
    
}

