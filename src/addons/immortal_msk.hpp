// Author:  Mario S. Könz <mskoenz@gmx.net>
// Date:    20.06.2013 16:55:46 EDT
// File:    immortal_msk.hpp

#ifndef __IMMORTAL_MSK_HEADER
#define __IMMORTAL_MSK_HEADER

#include "serialize/fserial.hpp"
#include <stdexcept>

namespace addon {
    class immortal_class {
        typedef size_t index_type;
    public:
        //------------------- ctors -------------------
        immortal_class(): file_switch_(2) //"invalid" at start
                        , index_(0) {
                            
            data_file_[0] = "immortal1.bin";
            data_file_[1] = "immortal2.bin";
            index_file_[0] = "immortal3.bin";
            index_file_[1] = "immortal4.bin";
            d_arch_[0].set_name(data_file_[0]);
            d_arch_[1].set_name(data_file_[1]);
        }
        //------------------- info -------------------
        bool available() {
            bool d1ex = exists(data_file_[0]);
            bool d2ex = exists(data_file_[1]);
            bool i1ex = exists(index_file_[0]);
            bool i2ex = exists(index_file_[1]);
            
            assert((i1ex and d1ex) or (!i1ex and !d1ex));
            assert((i2ex and d2ex) or (!i2ex and !d2ex));
            
            index_type i1 = i1ex ? try_to_read_index(0) : 0;
            index_type i2 = i2ex ? try_to_read_index(1) : 0;
            
            if(i1ex)
                if((i1 != 0 and i2 == 0) or (i1 >= i2)) {
                    file_switch_ = 0;
                    index_ = i1;
                    return true;
                }
            if(i2ex)
                if((i1 == 0 and i2 != 0) or (i2 >= i1)) {
                    file_switch_ = 1;
                    index_ = i2;
                    return true;
                }
            file_switch_ = 0;
            return false;
        }
        index_type try_to_read_index(uint8_t const & pos) {
            index_type res = 0;
            char end = 0;
            std::ifstream ifs(index_file_[pos], std::ios::in);
            ifs >> res;
            ifs >> end;
            if(end == 'e')
                return res;
            else
                return 0;
        }
        //------------------- ops -------------------
        void set_path(std::string const & path) {
            data_file_[0] = path + "/" + data_file_[0];
            data_file_[1] = path + "/" + data_file_[1];
            index_file_[0] = path + "/" + index_file_[0];
            index_file_[1] = path + "/" + index_file_[1];
            d_arch_[0].set_name(data_file_[0]);
            d_arch_[1].set_name(data_file_[1]);
        }
        template<typename T>
        void operator<<(T & t) {
            assert(file_switch_ < 2);
            d_arch_[file_switch_] << t;
        }
        template<typename T>
        void operator>>(T & t) {
            assert(file_switch_ < 2);
            d_arch_[file_switch_] >> t;
            file_switch_ = 1 - file_switch_;
        }
        index_type get_index(index_type const & defa = 0) {
            if(file_switch_ < 2) {
                return index_;
            }
            else {
                file_switch_ = 0;
                return defa;
            }
        }
        void write_next_index(index_type const & i) {
            index_ = i;
            std::ofstream ofs(index_file_[file_switch_]);
            ofs << index_;
            ofs << ' ';
            ofs << 'e';
            file_switch_ = 1 - file_switch_;
        }
        void reset() {
            remove(data_file_[0].c_str());
            remove(data_file_[1].c_str());
            remove(index_file_[0].c_str());
            remove(index_file_[1].c_str());
            file_switch_ = 0;
        }
    private:
        bool exists(std::string const & filename) {
            std::ifstream ifs(filename.c_str());
            return ifs.is_open();
        }
    private:
        uint8_t file_switch_;
        std::string data_file_[2];
        std::string index_file_[2];
        fserial_class d_arch_[2];
        
        index_type index_;
    } immortal;
}//end namespace addon

#endif //__IMMORTAL_MSK_HEADER
