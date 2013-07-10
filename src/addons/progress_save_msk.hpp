// Author:  Mario S. Könz <mskoenz@gmx.net>
// Date:    31.05.2013 15:07:05 EDT
// File:    progress_save_msk.hpp

#ifndef __PROGRESS_SAVE_MSK_HEADER
#define __PROGRESS_SAVE_MSK_HEADER

#include "color.hpp"

#include <map>
#include <fstream>
#include <sstream>
#include <utility>
#include <iostream>
#include <algorithm>
#include <stdexcept>
/*
//checkpoint.read(); //happens in ctor

for(unsigned j = checkpoint("j", j); j < 5; ++j)
    for(unsigned k = 0; k < 10; ++k)
        checkpoint.write();

checkpoint("j", j, 5) to set 5 as "0" value

//checkpoint.reset(); //happens in dtor

*/

namespace addon {
    class progress_save_class {
        typedef unsigned store_type;
        
    public:
        progress_save_class(): name_("checkpoint.txt") {
        }
        ~progress_save_class() {
            reset();
        }
        store_type operator()(std::string const & key, store_type & i, store_type defa = 0) {
            if(ref_dict_.find(key) != ref_dict_.end() and (ref_dict_[key]) != &i)
                throw std::runtime_error("progress_save_class::operator(): please use different keys for different references");
            else
                ref_dict_[key] = &i;
            
            //~ if(scope_ > scope_dict_[key]) {
                //~ 
            //~ }
            //~ else
                //~ 
            //~ ++scope_;
            
            if(dict_[key] > defa)
                std::swap(dict_[key], defa);
            
            return defa;
        }
        void write() {
            std::ofstream of_;
            //~ of_.open(name_.c_str(), std::ios_base::app);
            of_.open(name_.c_str());
            std::for_each(ref_dict_.begin(), ref_dict_.end(), 
                [&](std::pair<std::string const, store_type *> const & p) {
                    of_ << p.first << " " << *(p.second) << "   ";
                }
            );
            of_ << std::endl;
            of_.close();
        }
        void set_path(std::string const & path) {
            name_ = path + name_;
            read();
        }
        void reset() {
            remove(name_.c_str());
            dict_.clear();
        }
    private:
        void read() {
            std::string buffer = "";
            std::string last_line = "";
            std::ifstream if_;
            if_.open(name_.c_str());
            if(if_.is_open()) {
                do {
                    last_line = buffer;
                    getline(if_, buffer);
                } while(buffer.size() > 0 or !if_.eof());
                //~ std::cout << "ll " << last_line << std::endl;
                std::istringstream iss(last_line);
                
                std::cout << GREEN << "checkpoint file found" << NONE << std::endl;
                while(!iss.eof()) {
                    std::string key;
                    iss >> key;
                    if(key.size())
                        iss >> dict_[key];
                    else
                        break;
                    std::cout << GREEN << key << " = " << GREENB << dict_[key] << NONE << std::endl;
                }
                if_.close();
            }
            else
                std::cout << YELLOW << "no checkpoint file found" << NONE << std::endl;
        }
    private:
        std::string name_;
        std::map<std::string, store_type> dict_;
        
        static int scope_;
        std::map<std::string, int> scope_dict_;
        
        std::map<std::string, store_type *> ref_dict_;
    } checkpoint;
    
    int progress_save_class::scope_ = 0;
}//end namespace detail

#endif //__PROGRESS_SAVE_MSK_HEADER
