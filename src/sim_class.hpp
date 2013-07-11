// Author:  Mario S. Könz <mskoenz@gmx.net>
// Date:    18.04.2013 12:48:33 EDT
// File:    sim_class.hpp

#ifndef __SIM_CLASS_HEADER
#define __SIM_CLASS_HEADER

#include <jackknife.hpp>
#include <grid_class.hpp>
#include <shift_region_class.hpp>

#include <timer2_msk.hpp>
#include <random2_msk.hpp>
#include <accum_double.hpp>
#include <accum_simple.hpp>
#include <immortal_msk.hpp>
#include <bash_parameter3_msk.hpp>

#include <map>
#include <cmath>
#include <iostream>
#include <assert.h>
#include <algorithm>

//perimeter is documented in grid_class.hpp
namespace perimeter_rvb {
    ///  \brief holds update/measurement/checkpointsystem/rng in one place
    class sim_class {
        typedef typename grid_class::index_type index_type; ///< just forwarding from grid
        typedef typename grid_class::site_type site_type; ///< just forwarding from grid
        typedef addon::bash_parameter_class::map_type map_type; ///< just forwarding from bash_parameter
    public:
        ///  \brief the only constructor
        ///  
        ///  @param param is a map that contains the bash parameters
        ///  
        ///  param is used to specify all modifiable behavior of the simulation during runtime
        sim_class(map_type const & param):    param_(param)
                                            , H_(param_["H"])
                                            , L_(param_["L"])
                                            , grid_(H_, L_, std::vector<unsigned>(2, qmc::n_bonds == qmc::hex ? 2 : 0))
                                            , rngS_()
                                            , rngH_(H_)
                                            , rngL_(L_)
                                            {
                                                
            shift_region_class sr_(param_["shift"]);
            grid_.set_shift_region(sr_);
            
            //uncomment for negative "vortex" init in the triangular case (fails for sqr and hex)
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
        ///  \brief just forwards the two bond update to the grid with a random tile (tri only)
        bool two_bond_update(index_type i, index_type j, state_type state) {
            if(qmc::n_bonds == qmc::tri)
                return grid_.two_bond_update_intern(i, j, state, int(rngS_() * 3));
            else
                return grid_.two_bond_update_intern(i, j, state, 0);
        }
        ///  \brief just forwards the two bond update to the grid
        ///  
        ///  @param t specifies the tile. With this version one can manually update tiles (tri only)
        ///  
        ///  Two two_bond_update fcts are better bc there would be an if if we did only one version with a default for t
        bool two_bond_update(index_type i, index_type j, state_type state, int t) {
            if(qmc::n_bonds == qmc::tri)
                return grid_.two_bond_update_intern(i, j, state, t);
            else
                return grid_.two_bond_update_intern(i, j, state, 0);
        }
        ///  \brief changes the spin of the loops
        ///  
        ///  Decides at random (50:50) for every loop if all spins in the loop should be flipped or not
        void spin_update() {
            grid_.init_loops();
            
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
                                grid_.follow_loop_tpl(&s, bra, 
                                    [&](site_type * next) {
                                        next->check[bra] = true;
                                        next->spin[bra] = qmc::invert_spin - next->spin[bra];
                                    }
                                );
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
        ///  \brief updates bonds and spins
        ///  
        ///  does H*L update atempts on random tiles for each state followed by a spin_update
        void update() {
            grid_.set_shift_mode(qmc::no_shift);
            
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
            grid_.copy_to_ket(); //bc spins have changed
            grid_.clear_tile_spin(); //bc spins have changed
        }
        ///  \brief measures wanted properties
        ///  
        ///  just add your own data_["your_tag"] << your_value; to measure something
        void measure() {
            //=================== preswap zone ===================
            double loops = grid_.n_loops();
            
            data_["loops"] << loops;
            data_["overlap"] << pow(2.0, loops - 2*H_*L_* .5 );
            //=================== swap zone ===================
            grid_.set_shift_mode(qmc::ket_swap);
            grid_.init_loops();
            //grid_.eco_init_loops(); one could do a more efficient init_loop
            
            data_["sign"] << (grid_.sign() == 1);
            data_["neg_loops"] << grid_.n_neg_loops() / (double)grid_.n_loops();
            data_["swap_loops"] << grid_.n_loops();
            data_["swap_overlap"] << pow(2.0, int(grid_.n_loops()) - loops);
            data_["mean_for_error"] << pow(2.0, int(grid_.n_loops()) - loops);
            
            //=================== back to preswap ===================
            grid_.set_shift_mode(qmc::ket_preswap);
        }
        ///  \brief writes the the bins into the mean_file
        ///  
        ///  @param bins is a vector with values that need to go into the mean_file
        ///  @param mean_file is the file where the values are appended to
        ///  
        ///  since writing files is slow, we want to collect some values before writing. That why there's a vector
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
        ///  \brief core of the simulation
        ///  
        ///  see in code comments for detailed explanation
        void run() {
            //------------------- init state -------------------
            grid_.set_shift_mode(qmc::ket_preswap);
            //------------------- init timer -------------------
            addon::timer_class<addon::data> timer(param_["term"] + param_["sim"], param_["res"]);
            //set the labels for the later write
            timer.set_names("seed"
                          , "H"
                          , "sign"
                          , "sim"
                          , "x"
                          , "preswap_entropy"
                          , "neg_loops"
                          , "loop_time[us]"
                          , "entropy"
                          , "error"
                          ); //can take maximally 10 arguments
            timer.set_comment("measurement"); //optional, only shows in print not write
            
            //where the mean bins are
            std::string mean_file = (std::string)param_["prog_dir"] + "/mean.txt";
            
            std::vector<double> bins;
            
            //if -fix is found in the bash arguments it will just recalculate the jackknife
            if(param_.find("fix") == param_.end()) {
                //if -del is found in the bash arguments all progress/mean files will be deleted
                if(param_.find("del") != param_.end()) {
                    addon::immortal.reset();
                    remove(mean_file.c_str());
                    remove((std::string(param_["prog_dir"]) + "/state.txt").c_str()); //timer...
                }
                if(addon::immortal.available()) { //else if there are progress files around they get loaded
                    std::cout << GREENB << "load data at index " << addon::immortal.get_index() << NONE << std::endl;
                    addon::immortal >> (*this);
                }
                else {//otherwise the thermalization begins normally
                    //------------------- therm -------------------
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
                    
                    if((i & ((1lu<<10) - 1)) == ((1lu<<10) - 1)) { //all 1024 one mean gets added to the bins vector
                        
                        bins.push_back(data_["mean_for_error"].mean());
                        data_["mean_for_error"] = accumulator_double();
                        
                        if((i & ((1lu<<14) - 1)) == ((1lu<<14) - 1)) { //all 16*1024 the progress is saved and the bins written
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
            
            
            auto jack = jackknife(mean_file); //calc jackknife from mean_file
            
            timer.write(addon::global_seed.get()
                        , H_
                        , data_["sign"].mean()
                        , param_["sim"]
                        , param_["g"]
                        , -std::log(data_["swap_overlap"].mean())
                        , data_["neg_loops"].mean()
                        , timer.loop_time()
                        , jack.first
                        , jack.second
                        ); //can take maximally 10 arguments
        }
        ///  \brief just printing the data in the accumulators
        void present_data() {
            std::for_each(data_.begin(), data_.end(), 
                [&](std::pair<std::string const, accumulator_double> & p) {
                    std::cout << p.first << ": " << p.second << std::endl;
                }
            );
            std::cout << "S2 = " << -std::log(data_["swap_overlap"].mean()) << std::endl;
            std::cout << "accept " << int(accept_.mean() * 100) << "%" << std::endl;
        }
        ///  \brief small helper for simuviz_frame
        int simuviz_get_bond(site_type const & s, bond_type const & dir, bool const & spin_up) {
            bond_type bra_bond = s.bond[0];
            bond_type ket_bond = s.bond[qmc::invert_state - 0];
            
            if(bra_bond == dir or ket_bond == dir)
                return spin_up; //1==color during bond_updates / 0==color during spin_update
            else
                return 2; //2==invisible
        }
        ///  \brief writes a frame for a simuviz visualization
        ///  
        ///  @param spin_up says if the frame is taken during a spin_update or bond_update
        ///  
        ///  Frames are only taken, if the #define SIMUVIZ_FRAMES 50 is defined before including this header.
        ///  SIMUVIZ_FRAMES specifies how many frames are recorded. One needs to change the hex/sqr/tri_info.txt file
        ///  and adapt the size and the amount of frame (first two lines). The destination is whatever param_["vis"] returns
        void simuviz_frame(bool const & spin_up = false) {
            
            static int frame = 0;
            ++frame;
            if(frame > 50)
                return;
            DEBUG_VAR(frame)
            
            std::ofstream ofs;
            ofs.open(std::string(param_["vis"]), std::ios::app);
            
            ofs << std::endl << std::endl;
            
            if(qmc::n_bonds == qmc::tri) {
                for(unsigned i = 0; i < H_; ++i) {
                    for(unsigned j = 0; j < L_; ++j) {
                        ofs << grid_(i, j).spin[0] << " ";
                        ofs << simuviz_get_bond(grid_(i, j), qmc::up, spin_up) << " ";
                        ofs << simuviz_get_bond(grid_(i, j), qmc::diag_up, spin_up) << " ";
                        ofs << simuviz_get_bond(grid_(i, j), qmc::right, spin_up) << "  ";
                    }
                    ofs << std::endl;
                }
            } else if(qmc::n_bonds == qmc::hex) {
                for(unsigned i = 0; i < H_; ++i) {
                    for(unsigned j = 0; j < L_; j+=2) {
                        if(i % 2 == 0) {
                            ofs << grid_(i, j).spin[0] << " 2 2 " 
                               << simuviz_get_bond(grid_(i, j), qmc::up, spin_up) << " 2 2 2 2  " 
                               << grid_(i, j+1).spin[0] << " 2 " 
                               << simuviz_get_bond(grid_(i, j+1), qmc::up, spin_up) << " 2 " 
                               << simuviz_get_bond(grid_(i, j+1), qmc::hori, spin_up) << " 2 2 2  2 2 2 2 2 2 2 2  ";
                        }
                        else {
                            ofs << "2 " 
                               << grid_(i, j).spin[0] << " 2 2 2 " 
                               << simuviz_get_bond(grid_(i, j), qmc::up, spin_up) << " 2 " 
                               << simuviz_get_bond(grid_(i, j), qmc::hori, spin_up) << "  2 2 2 2 2 2 2 2  2 " 
                               << grid_(i, j+1).spin[0] << " 2 2 2 2 " 
                               << simuviz_get_bond(grid_(i, j+1), qmc::up, spin_up) << " 2  ";
                        }
                    }
                    ofs << std::endl;
                }
            }
            
            ofs.close();
        }
        ///  \brief just returns a reference to the grid
        grid_class & grid() {
            return grid_;
        }
        ///  \brief for the checkpoints
        ///  
        ///  this function is used by the serializer to get and set this object
        template<typename Archive>
        void serialize(Archive & ar) {
            ar & rngS_;
            ar & rngH_;
            ar & rngL_;
            ar & accept_;
            ar & grid_;
            ar & data_;
        }
    private:
        map_type param_;    ///< the parameter with all the settings
        const unsigned H_;  ///< height
        const unsigned L_;  ///< length
        grid_class grid_;   ///< the actual grid
        addon::random_class<double, addon::mersenne> rngS_; ///< spin/tile-random source
        addon::random_class<int, addon::mersenne> rngH_;    ///< H-random source
        addon::random_class<int, addon::mersenne> rngL_;    ///< L-random source
        
        std::map<std::string, accumulator_double> data_;    ///< all measurements are stored in here
        accumulator_simple accept_; ///< measures the update acceptance for bond_updates
    };
}
#endif //__SIM_CLASS_HEADER
