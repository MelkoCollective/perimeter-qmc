// Author:  Mario S. Könz <mskoenz@gmx.net>
// Date:    17.04.2013 10:45:13 EDT
// File:    grid_class.hpp

#ifndef __GRID_CLASS_HEADER
#define __GRID_CLASS_HEADER

#include <site_struct.hpp>

#include <boost/integer.hpp>
#include <boost/multi_array.hpp>

#include <map>
#include <memory>
#include <vector>
#include <iomanip>
#include <iostream>
#include <algorithm>

///  \brief where all the physics happens
///  
///  
namespace perimeter_rvb {
    ///  \brief grid class with site_struct as sites
    ///
    ///  via the defines in conf.hpp one can specified, what grid-type is chosen 
    ///  and how many states the grid has
    class grid_class { //grid is shorter than lattice ; 
    public:
        ///  the type of the sites
        typedef site_struct site_type;
    private:
        ///  the array where the sites are stored in
        template<typename U> 
        using array_type = boost::multi_array<U, 2>;
        
        ///  the used vector_type
        template<typename U> 
        using vector_type = std::vector<U>;
    public:
        ///  normally a size_t
        typedef typename vector_type<site_type>::size_type index_type;
        
        ///  \brief the only constructor
        ///  
        ///  @param H is the height of the grid
        ///  @param L is the length of the grid
        ///  @param init is a vector that can be as long as n_bra (shorter is fine too) that specifies the initialisation, has a default
        ///  
        ///  For sqr: init 2 is illegal
        ///  
        ///  For hex: init 2 is the only one that results in an updateable config
        ///  
        ///  bra and ket states are identical at the beginning
        grid_class(unsigned const H, unsigned const L, std::vector<unsigned> init = std::vector<unsigned>(qmc::n_bra, 0)): 
                alternator_(qmc::start_state)
              , H_(H)
              , L_(L)
              , grid_(boost::extents[H_][L_])
              , n_loops_(0)
              , shift_mode_(qmc::no_shift) {
            
            //just make sure that the input is sensible
            assert(H_>0);
            assert(L_>0);
            assert(qmc::n_states > 0);
            assert(H_%2 == 0);
            assert(L_%2 == 0);
            assert(qmc::n_states%2 == 0);
            assert(qmc::n_bonds == qmc::tri or qmc::n_bonds == qmc::sqr or qmc::n_bonds == qmc::hex);
            
            if(qmc::n_bonds == qmc::hex) {
                if(H_ % 3 != 0 or L_%3 != 0)
                    throw std::runtime_error("L and H must be divisible by 3 and 2 for the hex grid");
            }
            
            //if init is too short it gets filled up with zeros
            while(init.size() < qmc::n_bra)
                init.push_back(0);
            
            init_grid(init);
            
            init_tile();
            
            init_loops();
            
            copy_to_ket();
        }
        ///  \brief resets the check flag for all sites in all states
        ///  
        ///  Whenever an operation is performed that concernes the whole grid, visited sites will be flaged as checked.
        ///  At the end of the operation (e.g. init_loops) one has to clear this flags
        void clear_check(){
            std::for_each(begin(), end(), 
                [&](site_type & s) {
                    s.check = check_type();
                }
            );
        }
        ///  \brief 
        ///  
        ///  @param start marks the site from where the loop should be followed
        ///  @param bra stands for the "layer" (like a z-coordinate). Shouldn't be changed by this function is all works well (you return to the same point if you followed the loop)
        ///  @param fct is called for each site on the loop. fct(site_type *) must be legal in order to work.
        ///  
        ///  just traverses the loop where start is in an calls fct for all site it visites
        template<typename F>
        void follow_loop_tpl(site_type * const start, state_type & bra, F fct) {
            state_type old_bra = bra;
            site_type * next = start;
            
            do {
                fct(next);
                next = next_in_loop(next, bra);
            } while(next != start or bra != old_bra);
        }
        
        ///  \brief just forwards the work to tile_update
        ///  
        ///  @param i is the height coordinate
        ///  @param j is the length coordinate
        ///  @param state specifies in what bra or ket the update should be tried
        ///  @param tile (only relevant for triangular) specifies the tile
        ///  
        ///  The function returns true if the update was successful, false otherwise
        bool two_bond_update_intern(unsigned const & i, unsigned const & j, state_type const & state, unsigned const & tile) {
            return grid_[i][j].tile_update(state, tile);
        }
        ///  \brief reset the spin-checked-flags on the tiles
        ///  
        ///  After a spinupdate, the tiles have to be checked again if they are now or still updateable.
        ///  Clearing this bit will not check for that, but enable the check if a tile is visited
        void clear_tile_spin() {
            for(state_type state = qmc::start_state; state < qmc::n_states; ++state) {
                std::for_each(begin(), end(), 
                    [&](site_type & s) {
                        for(unsigned i = 0; i < tile_type::tile_per_site; ++i) {
                            CLEAR_BIT(s.tile[state][i].alpha, qmc::clear)
                        }
                    }
                );
            }
        }
        ///  \brief change from no swap to preswap or swap
        ///  
        ///  @param new_mode is the new desired shift_mode
        ///  
        ///  This function just sets the internal shift_mode_
        void set_shift_mode(shift_type const & new_mode) {
            site_type::shift_mode_print = new_mode; //just for printing
            shift_mode_ = new_mode;
        }
        ///  \brief writes the "shift matrix" to the sites
        ///  
        ///  Any object that has an operator()(size_t, size_t, size_t) and is equal or larger that the
        ///  works in here. The information is just transfered to the sites
        template<typename T> //T must suport an operator()(size_t, size_t, size_t)
        void set_shift_region(T const & region) {
            for(shift_type shift_mode = qmc::start_shift; shift_mode != qmc::n_shifts; ++shift_mode)
                for(index_type i = 0; i < H_; ++i)
                    for(index_type j = 0; j < L_; ++j)
                        grid_[i][j].shift_region[shift_mode] = region(shift_mode, i, j);
        }
        ///  \brief copies spins from bra to ket
        ///  
        ///  In case of no shift, the ket_spins are the same as the corresponding bra_spins.
        ///  If a shift is set, it will permute the spins accordingly. It's just a trick
        ///  for a nicer implementation, so that the two_bond_update doesn't have to "jump"
        ///  between layers (states)
        void copy_to_ket() {
            if(shift_mode_ == qmc::no_shift) {
                for(state_type bra = qmc::start_state; bra < qmc::n_bra; ++bra) {
                    std::for_each(begin(), end(), 
                        [&](site_type & s) {
                            state_type const ket = qmc::invert_state - bra;
                            s.spin[ket] = s.spin[bra];
                        }
                    );
                }
            }
            else {
                for(state_type bra = qmc::start_state; bra < qmc::n_bra; ++bra) {
                    std::for_each(begin(), end(), 
                        [&](site_type & s) {
                            state_type const ket = qmc::invert_state - bra;
                            state_type effective_bra = bra + (qmc::n_bra - s.shift_region[shift_mode_]);
                            if(effective_bra >= qmc::n_bra)//lazy boundary for now
                                effective_bra -= qmc::n_bra;
                            s.spin[ket] = s.spin[effective_bra];
                        }
                    );
                }
            }
        }
        
        ///  \brief return site at position (i / j)
        site_type & operator()(index_type const i, index_type const j) {
            return grid_[i][j];
        }
        ///  \brief return site at position (i / j) const version
        site_type const & operator()(index_type const i, index_type const j) const {
            return grid_[i][j];
        }
        ///  \brief initializes the loop structure and counts them
        ///  
        ///  the sign of the total config as well as the number of "negative" loops is also captured
        void init_loops() {
            n_loops_ = 0;
            n_neg_loops_ = 0;
            sign_ = +1;
            int subsign = +1;
            
            for(state_type bra = qmc::start_state; bra < qmc::n_bra; ++bra) //all transition graphs
                std::for_each(begin(), end(), // all sites in a transition graph
                    [&](site_type & s) {
                        if(s.check[bra] == false) { //only if not already visited
                            subsign = +1;
                            alternator_ = bra; //must be bra, not ket, see subsign *= -1 below
                            auto old_bra = bra;
                            site_type::last_dir = 0;
                            
                            follow_loop_tpl(&s, bra, 
                                [&](site_type * next){
                                    next->check[bra] = true;
                                    next->loop[bra] = n_loops_;
                                    assert(alternator_ == bra or alternator_ == qmc::invert_state - bra);
                                    
                                    if(alternator_ == bra) {
                                        if(site_type::last_dir >= qmc::middle) //will never happen at the beginning since last_dir == 0
                                            subsign *= -1;
                                    }
                                    else {
                                        if(site_type::last_dir < qmc::middle)
                                            subsign *= -1;
                                    }
                                }
                            );
                            if(site_type::last_dir >= qmc::middle) //here we need the bra case
                                subsign *= -1;
                            
                            assert(old_bra == bra); //we must end in the same layer we started
                            
                            if(subsign == -1) {
                                ++n_neg_loops_;
                                sign_ *= -1;
                            }
                            
                            ++n_loops_;
                        }
                    }
                );
            clear_check();
        }
        ///  \brief just returns internal n_loops_
        loop_type const & n_loops() const {
            return n_loops_;
        }
        ///  \brief just returns internal n_neg_loops_
        loop_type const & n_neg_loops() const {
            return n_neg_loops_;
        }
        ///  \brief just returns internal sign_
        int const & sign() const {
            return sign_;
        }
        //=================== print and iterate ===================
        ///  \brief prints the sites accordung to the site.print fct
        ///  
        ///  @param state is vector that specifies what states should be printed
        ///  @param os is the print destination
        void print(std::vector<state_type> state = std::vector<state_type>(), std::ostream & os = std::cout) const {
            for(index_type i = 0; i < state.size(); ++i) {
                state_type bra = state[i];
                os << "state nr: " << bra << std::endl;
                for(index_type i = 0; i < H_; ++i) {
                    for(index_type j = 0; j < L_; ++j) {
                        os << std::setw(3);
                        grid_[i][j].print(bra, os);
                        os << " ";
                    }
                    os << std::endl;
                }
            }
        }
        ///  \brief prints the sites accordung to the site.string_print fct
        ///  
        ///  @param state is vector that specifies what states should be printed
        ///  @param flag specifies how the states should be printed
        ///  @param os is the print destination
        ///  
        ///  If states contains a 1 (i.e. first bit is set) the transition graph is printed
        ///  
        ///  If states contains a 2 (i.e. second bit is set) the bra state is printed
        ///  
        ///  If states contains a 4 (i.e. third bit is set) the corresponding ket state is printed
        void print_all(std::vector<state_type> state = std::vector<state_type>(), unsigned flags = 1, std::ostream & os = std::cout) const {
            if(state.size() == 0) {
                for(state_type bra = qmc::start_state; bra != qmc::n_bra; ++bra) {
                    state.push_back(bra);
                }
            }
            for(index_type i = 0; i < state.size(); ++i) {
                state_type bra = state[i];
                os << "state nr: " << bra << std::endl;
                const unsigned kmax = 3;
                
                array_type<std::string> s(boost::extents[kmax * H_][5]);
                vector_type<std::string> in;
                if((flags & 1) == 1) { //print tragsition graph
                    site_struct::print_alternate = 1;
                    for(index_type i = 0; i < H_; ++i) {
                        for(index_type j = 0; j < L_; ++j) {
                            in = grid_[i][j].string_print(L_, bra, 0);
                            for(index_type k = 0; k < kmax; ++k) {
                                s[i * kmax + k][0] += in[k];
                                s[i * kmax + k][1] = "      ";
                            }
                        }
                    }
                }
                
                if((flags & 2) == 2) { //print bra
                    site_struct::print_alternate = 1;
                    for(index_type i = 0; i < H_; ++i) {
                        for(index_type j = 0; j < L_; ++j) {
                            in = grid_[i][j].string_print(L_, bra, 1);
                            for(index_type k = 0; k < kmax; ++k) {
                                s[i * kmax + k][2] += in[k];
                                s[i * kmax + k][3] = "      ";
                            }
                        }
                    }
                }
                if((flags & 4) == 4) { //print ket
                    site_struct::print_alternate = 1;
                    for(index_type i = 0; i < H_; ++i) {
                        for(index_type j = 0; j < L_; ++j) {
                            in = grid_[i][j].string_print(L_, bra, 2);
                            for(index_type k = 0; k < kmax; ++k) {
                                s[i * kmax + k][4] += in[k];
                            }
                        }
                    }
                }
                for(index_type i = 0; i < H_ * kmax; ++i) {
                    if(qmc::n_bonds == qmc::tri)
                        for(index_type j = 0; j < H_*kmax - i; ++j)
                            os << " ";
                        
                    for(index_type j = 0; j < 5; ++j)
                        os << s[i][j];
                    os << std::endl;
                }
            }
        }
        ///  \brief sub_iterator
        ///  
        ///  can be used if one wants to traverse the multi_array with just one loop and not 2 loops (2D)
        site_type * begin() {
            return grid_.data();
        }
        ///  \brief sub_iterator
        ///  
        ///  marks the end of the multi_array elements
        site_type * end() {
            return grid_.data() + grid_.num_elements();
        }
        ///  \brief for the checkpoints
        ///  
        ///  this function is used by the serializer to get and set this object
        template<typename Archive>
        void serialize(Archive & ar) {
            ar & n_loops_;
            ar & alternator_;
            ar & shift_mode_;
            ar & grid_;
        }
    private:
        ///  \brief initializes the periodic neighbor structur as well as the initial state
        ///  
        ///  @param init is a vector that is as long as n_bra that specifies the initialisation
        ///  
        ///  is only used by the constructor thus private
        void init_grid(std::vector<unsigned> const init) {
            int state = 0;
            
            for(state_type bra = qmc::start_state; bra < qmc::n_bra; ++bra) {
                state = 0;
                std::for_each(begin(), end(), 
                    [&](site_type & s) {
                        state_type ket = qmc::invert_state - bra;
                        s.spin[bra] = (state + state / L_)%2 == 0 ? qmc::beta : qmc::alpha;
                        s.spin[ket] = s.spin[bra];
                        
                        if(init[bra] == 0) {
                            if(qmc::n_bonds == qmc::hex) {
                                s.bond[bra] = qmc::hori;
                                s.bond[ket] = qmc::hori;
                            }
                            else {
                                s.bond[bra] = (state%2==0 ? qmc::right:qmc::left);
                                s.bond[ket] = (state%2==0 ? qmc::right:qmc::left);
                            }
                        }
                        else if(init[bra] == 1) {
                            s.bond[bra] = (state/L_%2==0 ? qmc::down:qmc::up);
                            s.bond[ket] = (state/L_%2==0 ? qmc::down:qmc::up);
                        }
                        else if(init[bra] == 2) {
                            if(qmc::n_bonds == qmc::hex) {
                                s.bond[bra] = (state/L_%3==0 ? qmc::down: (state/L_%3==2 ? qmc::hori : qmc::up));
                                s.bond[ket] = (state/L_%3==0 ? qmc::down: (state/L_%3==2 ? qmc::hori : qmc::up));
                            }
                            else if(qmc::n_bonds == qmc::tri) {
                                s.bond[bra] = (state/L_%2==0 ? qmc::diag_down:qmc::diag_up);
                                s.bond[ket] = (state/L_%2==0 ? qmc::diag_down:qmc::diag_up);
                            }
                        }
                        s.loop[bra] = 1-(state + state / L_)%2; //important for tile_init hex
                    ++state;
                    }
                );
            }
            //initialising the neighbor structure
            for(unsigned i = 0; i < H_; ++i) {
                for(unsigned j = 0; j < L_; ++j) {
                    grid_[i][j].neighbor[qmc::up] = &grid_[(i+H_-1)%H_][j];
                    grid_[i][j].neighbor[qmc::down] = &grid_[(i+1)%H_][j];
                    grid_[i][j].neighbor[qmc::me] = &grid_[i][j]; //me pointer points to opject itself (for monomers)
                    
                    //periodic boundaries happen here
                    if(qmc::n_bonds == qmc::hex) {
                        grid_[i][j].neighbor[qmc::hori] = &grid_[i][(j+L_ + 1 - 2*((i+j)%2) )%L_];
                    }
                    else {
                        grid_[i][j].neighbor[qmc::left] = &grid_[i][(j+L_-1)%L_];
                        grid_[i][j].neighbor[qmc::right] = &grid_[i][(j+1)%L_];
                    }
                    
                    if(qmc::n_bonds == qmc::tri) {
                        grid_[i][j].neighbor[qmc::diag_up] = &grid_[(i+H_-1)%H_][(j+L_-1)%L_];
                        grid_[i][j].neighbor[qmc::diag_down] = &grid_[(i+1)%H_][(j+1)%L_];
                    }
                }
            }
        }
        ///  \brief initializes the tile(s) for each site
        ///  
        ///  just calls set_info for all tiles
        void init_tile() {
            for(state_type state = qmc::start_state; state < qmc::n_states; ++state) {
                std::for_each(begin(), end(), 
                    [&](site_type & s) {
                        for(unsigned i = 0; i < tile_type::tile_per_site; ++i)
                            s.tile[state][i].set_info(&s, state, i);
                    }
                );
            }
        }
        ///  \brief during the loop update the next site in the loop is returned by this fct
        ///  
        ///  @param in is the entering site for which the neighbor is searched
        ///  @param bra is the current layer (imagine it like a z-coordinate). It can be changed by this fct
        ///  
        ///  Changes the alternator_ and bra if a "jump" occures
        site_type * next_in_loop(site_type * const in, state_type & bra) {
            alternator_ = qmc::invert_state - alternator_;
            return in->loop_partner(alternator_, bra, shift_mode_); //alternator can be changed, as well as bra
        }
    public:
        state_type alternator_; ///< alternates between bra and ket to "create" the transition graph
    private:
        unsigned const H_;           ///<height
        unsigned const L_;           ///<length
        array_type<site_type> grid_; ///< the actual grid
        
        loop_type n_loops_;     ///< amount of loops in the transition graph
        loop_type n_neg_loops_; ///< amount of "negative" loops in the transition graph
        int sign_;              ///< sign of complete config (-1 is n_neg_loops_ is odd / +1 else)
        shift_type shift_mode_; ///< the current shift mode (no swap, preswap, swap)
    };
}//end namespace perimeter_rvb
#endif //__GRID_CLASS_HEADER
