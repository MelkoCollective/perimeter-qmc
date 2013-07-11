// Author:  Mario S. Könz <mskoenz@gmx.net>
// Date:    14.06.2013 11:13:51 EDT
// File:    jackknife.hpp

#ifndef __JACKKNIFE_HEADER
#define __JACKKNIFE_HEADER

#include <vector>
#include <fstream>
#include <assert.h>
#include <algorithm>

namespace perimeter_rvb {
    std::pair<double, double> jackknife(std::string name) {
        std::vector<double> data;
        std::vector<double> fjack;
        double d = 0;
        
        std::ifstream ifs;
        ifs.open(name);
        assert(ifs.is_open());
        while(ifs) {
            ifs >> d;
            data.push_back(d);
        }
        data.pop_back(); //bc of the last \n
        ifs.close();
        
        double sum = accumulate(data.begin(), data.end(), 0.0);
        double N_1 = data.size() - 1;
        
        for(unsigned i = 0; i < data.size(); ++i)
            fjack.push_back(-std::log((sum - data[i]) / N_1));
        
        sum = accumulate(fjack.begin(), fjack.end(), 0.0);
        double mean = sum / fjack.size();
        double sq_sum = std::inner_product(fjack.begin(), fjack.end(), fjack.begin(), 0.0);
        double stmdev = std::sqrt(N_1 * (sq_sum / fjack.size() - mean * mean));
        return std::make_pair(mean, stmdev);
    }
}//end namespace 
#endif //__JACKKNIFE_HEADER
