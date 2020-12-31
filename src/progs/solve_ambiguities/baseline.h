
/* 
 * File:   baseline.h
 * Author: corbin
 *
 * Created on 12. September 2017, 14:18
 */

#ifndef BASELINE_H
#define BASELINE_H

#include "matrix.h"
#include <vector>

struct Baseline
{
    std::vector<int> obs_idxs_in_ncfile;
    ivg::Matrix residuals;
    ivg::Matrix delays;
    ivg::Matrix delays_sigma;
    ivg::Matrix ambiguity_spacing;
    ivg::Matrix shift;
    ivg::Matrix integerAmbiguities;
    bool fixed;
};   


#endif /* BASELINE_H */

