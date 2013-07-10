// Author:  Mario S. Könz <mskoenz@gmx.net>
// Date:    20.06.2013 15:24:39 EDT
// File:    map_support.hpp

#ifndef __MAP_SUPPORT_HEADER
#define __MAP_SUPPORT_HEADER

#include <map>
#include <algorithm>

#include "../archive_enum.hpp"

namespace addon {
    template<typename Key, typename Val, typename Archive>
    void serialize(Archive & ar, std::map<Key, Val> & arg) {
        unsigned size = arg.size();
        ar & size;
        if(Archive::type == archive_enum::output) {
            std::for_each(arg.begin(), arg.end(), 
                [&](std::pair<Key, Val> p) {
                    ar & p.first;
                    ar & p.second;
                }
            );
        } else if(Archive::type == archive_enum::input) {
            arg.clear();
            Key key;
            Val val;
            for(unsigned i = 0; i < size; ++i) {
                ar & key;
                ar & val;
                arg[key] = val;
            }
        }
    }
}//end namespace addon
#endif //__MAP_SUPPORT_HEADER
