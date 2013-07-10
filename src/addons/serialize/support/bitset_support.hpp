// Author:  Mario S. Könz <mskoenz@gmx.net>
// Date:    20.06.2013 14:00:35 EDT
// File:    bitset_support.hpp

#ifndef __BITSET_SUPPORT_HEADER
#define __BITSET_SUPPORT_HEADER

#include <bitset>
#include <algorithm>

#include "../archive_enum.hpp"
#include "../byte_operation.hpp"

namespace addon {
    template<unsigned N, typename Archive>
    void serialize(Archive & ar, std::bitset<N> & arg) {
        unsigned idx = 0;
        uint64_t in = 0;
        if(Archive::type == archive_enum::input) {
            for(unsigned i = 0; i < N/64; ++i) {
                ar & in;
                for(uint8_t j = 0; j < 64; ++j) {
                    arg[idx] = util::read_bit(in, j);
                    ++idx;
                }
            }
            ar & in;
            for(uint8_t j = 0; idx < N; ++j) {
                arg[idx] = util::read_bit(in, j);
                    ++idx;
            }
        }
        else if(Archive::type == archive_enum::output) {
            for(unsigned i = 0; i < N/64; ++i) {
                for(uint8_t j = 0; j < 64; ++j) {
                    util::write_bit(in, j, arg[idx]);
                    ++idx;
                }
                ar & in;
            }
            in = 0;
            for(uint8_t j = 0; idx < N; ++j) {
                util::write_bit(in, j, arg[idx]);
                ++idx;
            }
            ar & in;
        }
    }
}//end namespace addon
#endif //__BITSET_SUPPORT_HEADER
