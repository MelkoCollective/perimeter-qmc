// Author:  Mario S. Könz <mskoenz@gmx.net>
// Date:    18.04.2013 12:48:33 EDT
// File:    sim_class.hpp

#ifndef __SIM_CLASS_HEADER
#define __SIM_CLASS_HEADER

#include <grid_class.hpp>
#include <shift_region_class.hpp>
#include <jackknife.hpp>

#include <random2_msk.hpp>
#include <timer2_msk.hpp>
#include <immortal_msk.hpp>
#include <bash_parameter3_msk.hpp>
#include <accum_double.hpp>
#include <accum_simple.hpp>

#include <map>
#include <algorithm>
#include <iostream>
#include <assert.h>
#include <cmath>

//perimeter is documented in grid_class.hpp
namespace perimeter_rvb {
    class sim_class {
        typedef typename grid_class::index_type index_type;
        typedef typename grid_class::site_type site_type;
        typedef addon::bash_parameter_class::map_type map_type;
    public:
        sim_class(map_type const & param):    param_(param)
                                            , H_(param_["H"])
                                            , L_(param_["L"])
                                            , grid_(H_, L_, std::vector<unsigned>(2, qmc::n_bonds == qmc::hex ? 2 : 0))
                                            , rngS_()
                                            , rngH_(H_)
                                            , rngL_(L_)
                                            {
                                                
            shift_region_class sr_(param_["shift"]);
            //~ sr_.set_grow(std::vector<bond_type>(1, qmc::right));
            //~ shift_region_class sr_(H_, L_, param_["spacing"]);
            //~ sr_.grow_partial(param_["g"]);
            //~ sr_.write(param_["shift"]);
            grid_.set_shift_region(sr_);
            
            //~ state_type state = 0;
            //~ if(qmc::n_bonds == qmc::tri) {
                //~ for(unsigned i = 0; i < H_; i += 4) {
                    //~ for(unsigned j = 0; j < L_; j+= 4) {
                        //~ two_bond_update(i+3, j+1, state, 1);
                        //~ two_bond_update(i+3, j+3, state, 1);
                        //~ two_bond_update(i+3, j+3, qmc::invert_state - state, 1);
                        //~ two_bond_update(i+3, j+2, state, 1);
                        //~ two_bond_update(i+2, j+1, qmc::invert_state - state, 1);
                        //~ grid_(i+2, j+1).spin[0] = qmc::invert_spin - grid_(i+2, j+1).spin[0];
                        //~ grid_(i+1, j+1).spin[0] = qmc::invert_spin - grid_(i+1, j+1).spin[0];
                        //~ grid_.clear_tile_spin();
                        //~ grid_.copy_to_ket();
                        //~ two_bond_update(i+2, j+1, qmc::invert_state - state, 0);
                        //~ grid_(i+2, j+2).spin[0] = qmc::invert_spin - grid_(i+2, j+2).spin[0];
                        //~ grid_(i+1, j+1).spin[0] = qmc::invert_spin - grid_(i+1, j+1).spin[0];
                        //~ grid_.clear_tile_spin();
                        //~ grid_.copy_to_ket();
                        //~ two_bond_update(i+1, j+1, state, 2);
                    //~ }
                //~ }
            //~ }
            
            grid_.clear_tile_spin();
            grid_.copy_to_ket();
        }
        
        bool two_bond_update(index_type i, index_type j, state_type state, int t = -1) {
            if(qmc::n_bonds == qmc::tri)
                if(t == -1)
                    return grid_.two_bond_update_intern(i, j, state, int(rngS_() * 3));
                else
                    return grid_.two_bond_update_intern(i, j, state, t);
            else
                return grid_.two_bond_update_intern(i, j, state, 0);
        }
        
        void spin_update() {
            
            grid_.init_loops();
            //~ DEBUG_MSG("start spin")
            for(state_type bra = qmc::start_state; bra < qmc::n_bra; ++bra) {
                grid_.alternator_ = bra;
                std::for_each(grid_.begin(), grid_.end(), 
                    [&](site_type & s) {
                        if(s.check[bra] == false)
                        {
                            if(rngS_() > .5) {
                                grid_.follow_loop_tpl(&s, bra, 
                                    [&](site_type * next) {
                                        next->check[bra] = true;
                                    }
                                );
                            }
                            else {
                                //~ std::cout << "follow change loop: " << s.loop[bra] << "   >";
                                grid_.follow_loop_tpl(&s, bra, 
                                    [&](site_type * next) {
                                        next->check[bra] = true;
                                        //~ std::cout << "."; std::cout.flush();
                                        next->spin[bra] = qmc::invert_spin - next->spin[bra];
                                    }
                                );
                                //~ std::cout << "done" << std::endl;
                                #ifdef SIMUVIZ_FRAMES
                                    simuviz_frame(1);
                                #endif //SIMUVIZ_FRAMES
                            }
                        }
                    }
                );
            }
            grid_.clear_check();
        }
        
        void update() {
            grid_.set_shift_mode(qmc::no_shift);
            
            grid_.clear_tile_spin();
            
            for(state_type state = qmc::start_state; state < qmc::n_states; ++state)
                for(unsigned i = 0; i < H_ * L_; ++i) {
                    bool ok = two_bond_update(rngH_(), rngL_(), state);
                    accept_ << ok;
                    #ifdef SIMUVIZ_FRAMES
                        if(ok)
                            simuviz_frame();
                    #endif //SIMUVIZ_FRAMES
                }
            
            grid_.set_shift_mode(qmc::ket_preswap);
            spin_update();
            grid_.copy_to_ket();
        }
        
        void measure() {
            //~ int sign = grid_.sign();
            double loops = grid_.n_loops();
            //~ double neg_loops = grid_.n_neg_loops();
            data_["loops"] << loops;
            data_["overlap"] << pow(2.0, loops - 2*H_*L_* .5 );
            grid_.set_shift_mode(qmc::ket_swap);
            grid_.init_loops();
            
            //~ if(sign != grid_.sign() and grid_.n_loops() > 12) {
                //~ std::cout << "\033[1;36m" << "---checkpoint 2---" << "\033[0m" << std::endl;
                //~ grid_.print_all({0, 1}, 7);
                //~ DEBUG_VAR(sign)
                //~ DEBUG_VAR(grid_.sign())
                //~ std::cin.get();
            //~ }
            
            //~ grid_.eco_init_loops();
            data_["sign"] << (grid_.sign() == 1);
            data_["neg_loops"] << grid_.n_neg_loops()/(double)grid_.n_loops();
            data_["swap_loops"] << grid_.n_loops();
            data_["swap_overlap"] << pow(2.0, int(grid_.n_loops()) - loops);
            data_["mean_for_error"] << pow(2.0, int(grid_.n_loops()) - loops);
            
            grid_.set_shift_mode(qmc::ket_preswap);
        }
        void write_bins(std::vector<double> & bins, std::string const & mean_file) {
            std::ofstream ofs;
            ofs.open(mean_file, std::ios::app);
            std::for_each(bins.begin(), bins.end(), 
                [&](double const & d) {
                    ofs << d << std::endl;
                }
            );
            ofs.close();
            bins.clear();
        }
        void run() {
            //------------------- init state -------------------
            grid_.set_shift_mode(qmc::ket_preswap);
            //------------------- init timer -------------------
            addon::timer_class<addon::data> timer(param_["term"] + param_["sim"], param_["res"]);
            timer.set_names("seed"
                          , "H"
                          , "sign"
                          , "sim"
                          , "x"
                          , "preswap_entropy"
                          //~ , "accept"
                          , "neg_loops"
                          , "loop_time[us]"
                          , "entropy"
                          , "error"
                          );
            timer.set_comment("test");
            
            std::string mean_file = (std::string)param_["prog_dir"] + "/mean.txt";
            
            std::vector<double> bins;
            
            if(param_.find("fix") == param_.end()) {
                
                if(param_.find("del") != param_.end()) {
                    addon::immortal.reset();
                    remove(mean_file.c_str());
                    remove((std::string(param_["prog_dir"]) + "/state.txt").c_str()); //timer...
                }
                if(addon::immortal.available()) {
                    std::cout << GREENB << "load data at index " << addon::immortal.get_index() << NONE << std::endl;
                    addon::immortal >> (*this);
                }
                else {
                    //------------------- term -------------------
                    std::cout << std::endl;
                    for(unsigned i = 0; i < param_["term"]; ++i) {
                        update();
                        timer.progress(i, param_["timer_dest"]);
                    }
                }
                //------------------- sim -------------------
                std::ofstream ofs;
                for(unsigned i = addon::immortal.get_index(0); i < param_["sim"]; ++i) {
                    
                    update();
                    measure();
                    timer.progress(param_["term"] + i, param_["timer_dest"]);
                    
                    if((i & ((1lu<<10) - 1)) == ((1lu<<10) - 1)) {
                        
                        bins.push_back(data_["mean_for_error"].mean());
                        data_["mean_for_error"] = accumulator_double();
                        
                        if((i & ((1lu<<14) - 1)) == ((1lu<<14) - 1)) {
                            //------------------- write out bins -------------------
                            write_bins(bins, mean_file);
                            //------------------- serialize config -------------------
                            addon::immortal << (*this);
                            addon::immortal.write_next_index(i + 1);
                        }
                    }
                }
                write_bins(bins, mean_file);
                timer.write_state(param_["term"] + param_["sim"]);
            }
            
            
            auto jack = jackknife(mean_file);
            
            timer.write(addon::global_seed.get()
                        , H_
                        , data_["sign"].mean()
                        , param_["sim"]
                        , param_["g"]
                        , -std::log(data_["swap_overlap"].mean())
                        //~ , accept_.mean()
                        , data_["neg_loops"].mean()
                        , timer.loop_time()
                        , jack.first
                        , jack.second
                        );
        }
        
        void present_data() {
            std::for_each(data_.begin(), data_.end(), 
                [&](std::pair<std::string const, accumulator_double> & p) {
                    std::cout << p.first << ": " << p.second << std::endl;
                }
            );
            std::cout << "S2 = " << -std::log(data_["swap_overlap"].mean()) << std::endl;
            std::cout << "accept " << int(accept_.mean() * 100) << "%" << std::endl;
        }
        
        int get_bond(site_type const & s, bond_type const & dir, bool const & spin_up) {
            bond_type bra = s.bond[0];
            bond_type ket = s.bond[qmc::invert_state - 0];
            
            if(bra == dir or ket == dir)
                return spin_up;
            else
                return 2;
        }
        void simuviz_frame(bool const & spin_up = false) {
            static int frame = 0;
            ++frame;
            if(frame > 50)
                return;
            DEBUG_VAR(frame)
            
            std::ofstream ofs;
            ofs.open(std::string(param_["prog_dir"]) + "/../../../SimuViz/Example/viz.txt", std::ios::app);
            
            ofs << std::endl << std::endl;
            
            if(qmc::n_bonds == qmc::tri) {
                for(unsigned i = 0; i < H_; ++i) {
                    for(unsigned j = 0; j < L_; ++j) {
                        ofs << grid_(i, j).spin[0] << " ";
                        ofs << get_bond(grid_(i, j), qmc::up, spin_up) << " ";
                        ofs << get_bond(grid_(i, j), qmc::diag_up, spin_up) << " ";
                        ofs << get_bond(grid_(i, j), qmc::right, spin_up) << "  ";
                    }
                    ofs << std::endl;
                }
            } else if(qmc::n_bonds == qmc::hex) {
                for(unsigned i = 0; i < H_; ++i) {
                    for(unsigned j = 0; j < L_; j+=2) {
                        if(i % 2 == 0) {
                            ofs << grid_(i, j).spin[0] << " 2 2 " 
                               << get_bond(grid_(i, j), qmc::up, spin_up) << " 2 2 2 2  " 
                               << grid_(i, j+1).spin[0] << " 2 " 
                               << get_bond(grid_(i, j+1), qmc::up, spin_up) << " 2 " 
                               << get_bond(grid_(i, j+1), qmc::hori, spin_up) << " 2 2 2  2 2 2 2 2 2 2 2  ";
                        }
                        else {
                            ofs << "2 " 
                               << grid_(i, j).spin[0] << " 2 2 2 " 
                               << get_bond(grid_(i, j), qmc::up, spin_up) << " 2 " 
                               << get_bond(grid_(i, j), qmc::hori, spin_up) << "  2 2 2 2 2 2 2 2  2 " 
                               << grid_(i, j+1).spin[0] << " 2 2 2 2 " 
                               << get_bond(grid_(i, j+1), qmc::up, spin_up) << " 2  ";
                        }
                    }
                    ofs << std::endl;
                }
            }
            
            ofs.close();
        }
        
        grid_class & grid() {
            return grid_;
        }
        template<typename Archive>
        void serialize(Archive & ar) {
            ar & rngS_;
            ar & rngH_;
            ar & rngL_;
            ar & accept_;
            ar & grid_;
            ar & data_;
        }
    //~ private:
        map_type param_;
        const unsigned H_;
        const unsigned L_;
        grid_class grid_;
        addon::random_class<double, addon::mersenne> rngS_;
        addon::random_class<int, addon::mersenne> rngH_;
        addon::random_class<int, addon::mersenne> rngL_;
        
        std::map<std::string, accumulator_double> data_;
        accumulator_simple accept_;
    };
}
#endif //__SIM_CLASS_HEADER
