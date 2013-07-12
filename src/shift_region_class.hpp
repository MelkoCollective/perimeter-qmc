// Author:  Mario S. Könz <mskoenz@gmx.net>
// Date:    06.05.2013 14:45:54 EDT
// File:    shift_region_class.hpp

#ifndef __SWAP_REGION_CLASS_HEADER
#define __SWAP_REGION_CLASS_HEADER

#include <string>
#include <fstream>
#include <iostream>

#include <site_struct.hpp>

namespace perimeter_rvb {
    ///  \brief just reads a shift file and behaves like a matrix
    class shift_region_class {
    public:
        ///  \brief only constructor
        ///  
        ///  @param filename is the input-file
        ///  
        ///  the format of the input file must follow this rules:
        ///  
        ///  (first the preswap, one and only one spaces between numbers are mandatory)
        ///  
        ///  0 0 
        ///  
        ///  0 0 
        ///  
        ///  (not even a space in this empty line!)
        ///  
        ///  1 0
        ///  
        ///  1 0
        ///  
        ///  (again no space in this empty line!)
        shift_region_class(std::string filename) {
            std::ifstream in(filename);
            std::string temp;
            
            if(in.is_open()) {
                unsigned consistent(0);
                
                stage1_.clear();
                stage1_.push_back(std::vector<std::string>());
                
                while(in) {
                    getline(in, temp);
                    if(temp.size() > 1) {
                        if(consistent == 0)
                            consistent = count_sites(temp);
                        if(count_sites(temp) == consistent) {
                            stage1_.back().push_back(temp);
                        }
                        else
                            throw std::runtime_error("inconsistent amount of sites in shift_region_class constructor");
                    }
                    else {
                        if(stage1_.back().size() > 0)
                            stage1_.push_back(std::vector<std::string>());
                    }
                };
                if(stage1_.back().size() == 0)
                    stage1_.pop_back();
                H_ = stage1_[0].size();
                L_ = consistent;
                N_ = stage1_.size();
            }
            else {
                std::stringstream ss;
                ss << "file: " << filename << " not found" << std::endl;
                throw std::runtime_error(ss.str());
            }
            convert_1_to_2();
        }
        ///  \brief prints the content
        ///  
        ///  @param flags determine what is printed
        ///  
        ///  if the flag contains 1 (first bit set) it will print the exact thing it read from the file.
        ///  If the flag contains 2 (second bit set) it will print what it parsed and what the matrix looks like
        void print(unsigned flags = 1) const {
            if((flags&1) == 1) {
                std::cout << "--------stage1-graphical--------" << std::endl;
                for(unsigned n = 0; n < stage1_.size(); ++n) {
                    std::cout << "state: " << n << std::endl;
                    for(unsigned i = 0; i < stage1_[n].size(); ++i) {
                        std::cout << "    " << stage1_[n][i] << std::endl;
                    }
                }
            }
            if((flags&2) == 2) {
                std::cout << "--------stage2-matrix--------" << std::endl;
                for(unsigned n = 0; n < N_; ++n) {
                    std::cout << "state: " << n << std::endl;
                    for(unsigned i = 0; i < H_; ++i) {
                        std::cout << "    ";
                        for(unsigned j = 0; j < L_; ++j) {
                            std::cout << (*this)(n, i, j) << " ";
                        }
                        std::cout << std::endl;
                    }
                }
            }
        }
        ///  \brief return matrix-elements with this operator
        unsigned operator()(unsigned const & n, unsigned const & i, unsigned const & j) const {
            return stage2_[n][i][j];
        }
    private:
        ///  \brief parses the input into a matrix of numbers
        void convert_1_to_2() {
            stage2_ = std::vector<std::vector<std::vector<unsigned>>>(N_, std::vector<std::vector<unsigned>>(H_, std::vector<unsigned>(L_, 0)));
            for(unsigned n = 0; n < N_; ++n) {
                for(unsigned i = 0; i < H_; ++i) {
                    for(unsigned j = 0; j < L_; ++j) {
                        stage2_[n][i][j] = (stage1_[n][i][2*j] - '0');
                    }
                }
            }
        }
        ///  \brief just makes sure that each line contains the same amount of sites (chars between '0'-'9')
        unsigned count_sites(std::string const & in) const {
            unsigned res(0);
            for(unsigned i = 0; i < in.size(); ++i) {
                if(in[i] >= '0' and in[i] <= '9')
                    ++res;
            }
            return res;
        }
    private:
        unsigned H_;    ///< height
        unsigned L_;    ///< length
        unsigned N_;    ///< number of shift regions (preswap/swap)
        std::vector<std::vector<std::string>> stage1_;  ///< the information that came from the file
        std::vector<std::vector<std::vector<unsigned>>> stage2_; ///< the parsed matrix
    };
}//end namespace perimeter_rvb

#endif //__SWAP_REGION_CLASS_HEADER
