// Author:  Mario S. Könz <mskoenz@gmx.net>
// Date:    17.04.2013 10:45:13 EDT
// File:    grid_class.hpp

#ifndef __GRID_CLASS_HEADER
#define __GRID_CLASS_HEADER

#include <site_struct.hpp>

#include <boost/integer.hpp>
#include <boost/multi_array.hpp>

#include <iostream>
#include <iomanip>
#include <vector>
#include <map>
#include <set>
#include <memory>
#include <algorithm>
#include <typeinfo>

///  \brief where all the physics happens
///  
///  
namespace perimeter {
    ///  \brief grid class with site_struct as sites
    ///
    ///  via the enums in site_struct.hpp can be specified, what grid-type is chosen 
    ///  and how many states pairs of states the grid has
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
        
        inline grid_class(unsigned const H, unsigned const L, std::vector<unsigned> init = std::vector<unsigned>(qmc::n_bra, 0)): 
                H_(H)
              , L_(L)
              , grid_(boost::extents[H_][L_])
              , n_loops_(0)
              , alternator_(qmc::start_state)
              , shift_mode_(qmc::no_shift) {
            assert(H_%2 == 0);
            assert(L_%2 == 0);
            assert(qmc::n_states%2 == 0);
            assert(H_>0);
            assert(L_>0);
            assert(qmc::n_states > 0);
            assert(qmc::n_bonds == qmc::tri or qmc::n_bonds == qmc::sqr or qmc::n_bonds == qmc::hex);
            
            if(qmc::n_bonds == qmc::hex) {
                if(H_ % 3 != 0 or L_%3 != 0)
                    throw std::runtime_error("L and H must be divisible by 3 for the hex grid");
            }
            
            while(init.size() < qmc::n_bra)
                init.push_back(0);
                
            init_grid(init);
            
            init_tile();
            
            init_loops();
            
            copy_to_ket();
        }
        void clear_check(){
            std::for_each(begin(), end(), 
                [&](site_type & s) {
                    s.check = check_type();
                }
            );
        }
        template<typename F>
        void follow_loop_tpl(site_type * const start, state_type & bra, F fct) {
            state_type old_bra = bra;
            site_type * next = start;
            
            do {
                fct(next);
                next = next_in_loop(next, bra);
            } while(next != start or bra != old_bra);
        }
        void two_bond_flip(site_type * const target, site_type * const old_partner, bond_type const & b, state_type const & state) {
            //target node shows in the dircetion of the neighbor with the same orientation
            target->bond[state] = b;
            //old partner does the same
            old_partner->bond[state] = b;
            //the new partner of the target bond shows in the targets direction
            target->neighbor[b]->bond[state] = qmc::invert_bond - b;
            //old partner of the new partner does the same
            old_partner->neighbor[b]->bond[state] = qmc::invert_bond - b;
        }
        //------------------- omega impl -------------------
        bool two_bond_update_intern_2(unsigned const & i, unsigned const & j, state_type const & state, unsigned const & tile) {
            return grid_[i][j].tile_update(state, tile);
        }
        bool two_bond_update_intern(unsigned const & i, unsigned const & j, state_type const & state) {
            //only for square
            assert(qmc::n_bonds == qmc::sqr);
            
            site_type & target = grid_[i][j];
            
            state_type bra = state;
            if(bra >= qmc::n_bra) //it's a ket
                bra = qmc::invert_state - state; //now it's a bra
            if(target.bond[state] < 3) { // it's perhaps updateable
                for(bond_type b = qmc::start_bond; b < 3; ++b) {
                    if(target.bond[state] == target.neighbor[b]->bond[state] and target.spin[state] == qmc::invert_spin - target.neighbor[b]->spin[state]) {
                        two_bond_flip(&target, target.partner(state), b, state);
                        //~ res = true;
                        //~ break;
                        return true;
                    }
                }
            }
            return false;
        }
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
        void set_shift_mode(shift_type const & new_mode) {
            site_type::shift_mode_print = new_mode;
            shift_mode_ = new_mode;
        }
        template<typename T> //T must suport an operator()(size_t, size_t, size_t)
        void set_shift_region(T const & region) {
            for(shift_type shift_mode = qmc::start_shift; shift_mode != qmc::n_shifts; ++shift_mode)
                for(index_type i = 0; i < H_; ++i)
                    for(index_type j = 0; j < L_; ++j)
                        grid_[i][j].shift_region[shift_mode] = region(shift_mode, i, j);
        }
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
        
        site_type & operator()(index_type const i, index_type const j) {
            return grid_[i][j];
        }
        site_type const & operator()(index_type const i, index_type const j) const {
            return grid_[i][j];
        }
        
        std::pair<bool, double> init_loops() { //returns true if sign == -1 and the amount of negative loops
            n_loops_ = 0;
            n_neg_loops_ = 0;
            sign_ = +1;
            int subsign = +1;
            
            std::stringstream ss;
            
            for(state_type bra = qmc::start_state; bra < qmc::n_bra; ++bra)
                std::for_each(begin(), end(), 
                    [&](site_type & s) {
                        if(s.check[bra] == false) {
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
                            
                            assert(old_bra == bra);
                            //~ std::cout << "loop " << n_loops_ << " done " << bra << std::endl;
                            if(subsign == -1) {
                                ++n_neg_loops_;
                                sign_ *= -1;
                            }
                            
                            ss.clear();
                            ++n_loops_;
                        }
                    }
                );
            clear_check();
            return std::make_pair(sign_ == -1, n_neg_loops_/n_loops_);
        }
        loop_type const & n_loops() const {
            return n_loops_;
        }
        loop_type const & n_neg_loops() const {
            return n_neg_loops_;
        }
        int const & sign() const {
            return sign_;
        }
        //=================== print and iterate ===================
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
        site_type * begin() {
            return grid_.data();
        }
        site_type * end() {
            return grid_.data() + grid_.num_elements();
        }
        template<typename Archive>
        void serialize(Archive & ar) {
            ar & n_loops_;
            ar & alternator_;
            ar & shift_mode_;
            ar & grid_;
        }
    private:
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
                    grid_[i][j].neighbor[qmc::me] = &grid_[i][j];
                    
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
        site_type * next_in_loop(site_type * const in, state_type & bra) {
            alternator_ = qmc::invert_state - alternator_;
            return in->loop_partner(alternator_, bra, shift_mode_); //alternator can be changed, as well as bra
        }
    private:
        unsigned const H_; ///<height
        unsigned const L_; ///<length
        array_type<site_type> grid_; ///< the actual grid
        
        loop_type n_loops_;
        loop_type n_neg_loops_;
        int sign_;
    public:
        state_type alternator_;
    private:
        shift_type shift_mode_;
    };
}//end namespace perimeter
#endif //__GRID_CLASS_HEADER
