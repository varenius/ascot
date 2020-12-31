/* 
 * File:   CableCor.h
 * Author: armin
 *
 * Created on 19. Oktober 2015, 09:19
 */

#ifndef CABLECOR_H
#define	CABLECOR_H

#include "date.h"

namespace ltn{

class CableCor {
public:
    CableCor();
    virtual ~CableCor();
    CableCor(ivg::Date date, double corr);
    
    // Getter Setter
    void SetCorr(double corr);
    double GetCorr() const;
    void SetDate(ivg::Date date);
    ivg::Date GetDate() const;
private:
    ivg::Date date;
    double corr;

};

}

#endif	/* CABLECOR_H */

