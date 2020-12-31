#ifndef PERFORM_ANALYSIS_H
#define	PERFORM_ANALYSIS_H

#include <QThread>

#include <iostream>
#include <cstdlib>
#include <libconfig.h++>

#include "session_inout.h"
#include "session.h"

enum threadtype {
    LOAD,
    INIT,
    RESIDUALS,
    EXPORT
};

class Ascotizer : public QThread
{
    Q_OBJECT
    
public:
    
    Ascotizer();
    
    Ascotizer(threadtype type, int session_idx, ivg::Session &session);
        
    int get_session_idx(){ return _idx; };
    
    threadtype get_type(){ return _type; };
    
    ivg::Session & get_session(){ return _session; };
    
    
protected:
    
    void run();
    
private:
    
    int _idx;
    threadtype _type;
    ivg::Session _session;
    
};

#endif	/* PERFORM_ANALYSIS_H */

