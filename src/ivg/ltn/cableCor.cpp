/* 
 * File:   CableCor.cpp
 * Author: armin
 * 
 * Created on 19. Oktober 2015, 09:19
 */

#include "cableCor.h"
#include "date.h"

namespace ltn{

CableCor::CableCor() {
}

CableCor::~CableCor() {
}

CableCor::CableCor(ivg::Date date, double corr){
    this->date = date;
    this->corr = corr;
}


// Getter Setter
void CableCor::SetCorr(double corr) {
    this->corr = corr;
}

double CableCor::GetCorr() const {
    return corr;
}

void CableCor::SetDate(ivg::Date date) {
    this->date = date;
}

ivg::Date CableCor::GetDate() const {
    return date;
}

}

