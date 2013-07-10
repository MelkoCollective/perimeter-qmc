// Author:  Mario S. Könz <mskoenz@gmx.net>
// Date:    20.06.2013 13:21:04 EDT
// File:    boost_multiarray_support.hpp

#ifndef __BOOST_MULTIARRAY_SUPPORT_HEADER
#define __BOOST_MULTIARRAY_SUPPORT_HEADER

#include <boost/integer.hpp>
#include <boost/multi_array.hpp>
#include <algorithm>

namespace addon {
    template<typename T, unsigned N, typename Archive>
    void serialize(Archive & ar, boost::multi_array<T, N> & arg) {
        for(T * i = arg.data(); i < arg.data() + arg.num_elements(); ++i)
            ar & *i;
    }
}//end namespace addon
#endif //__BOOST_MULTIARRAY_SUPPORT_HEADER
