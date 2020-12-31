#ifndef _EOPSERIES_H
#define	_EOPSERIES_H

#include "date.h"
#include "eop_series.h"
#include "ui_eopseries.h"

class EopSeries : public QDialog {
    
    Q_OBJECT
    
public:
    
    EopSeries();
    
    virtual ~EopSeries();
    
    void set_start_end_epochs(ivg::Date start, ivg::Date end);
    
    void set_eop_file(string name, string type, string path);
    
    ivg::Eop_series get_eop_series(string &name);
    
private:
    
    Ui::EopSeries _ui;
    
};

#endif	/* _EOPSERIES_H */
