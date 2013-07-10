// Author:  Mario S. Könz <mskoenz@gmx.net>
// Date:    20.06.2013 15:18:33 EDT
// File:    string_support.hpp

#ifndef __STRING_SUPPORT_HEADER
#define __STRING_SUPPORT_HEADER

#include <string>

#include "../archive_enum.hpp"

namespace addon {
    template<typename Archive>
    void serialize(Archive & ar, std::string & arg) {
        unsigned size = arg.size();
        ar & size;
        if(Archive::type == archive_enum::input)
            arg.resize(size);
        for(unsigned i = 0; i < size; ++i) {
            ar & arg[i];
        }
    }
}//end namespace addon
#endif //__STRING_SUPPORT_HEADER
