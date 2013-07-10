// Author:  Mario S. Könz <mskoenz@gmx.net>
// Date:    18.04.2013 11:08:01 EDT
// File:    site_struct.hpp

#ifndef __SITE_STRUCT_HEADER
#define __SITE_STRUCT_HEADER

#include <iostream>
#include <iomanip>
#include <sstream>
#include <vector>
#include <bitset>
#include <assert.h>

#include <boost/integer.hpp>
#include <addons/color.hpp>

#include <conf.hpp>

#define STRINGIFY_HELPER(x) #x
#define STRINGIFY(x) STRINGIFY_HELPER(x) //is needed, believe it or not ...

#define DEBUG_VAR(x) std::cout << "\033[1;31m" << "  DEBUG_VAR: " << "\033[0;31m" << #x << " = " << x << "\033[0m" << std::endl;
#define DEBUG_MSG(x) std::cout << "\033[1;31m" << "  DEBUG_MSG: " << "\033[0;31m" << x << "\033[0m" << std::endl;

#define SET_BIT(x, y) (x) |= (y);
#define CLEAR_BIT(x, y) (x) &= ~(y);

//perimeter is documented in grid_class.hpp
namespace perimeter {
    ///  \brief namespace for all compiletime information
    ///  
    ///  contains all enums that determine the grid type, the amount of states and also the spin-states
    namespace qmc {
        ///  \brief contains the information about the states
        ///  
        ///  n_states has to be even, since a transition graph always needs a bra and a ket.
        ///  the bra and the corresponding ket have to be arranged in such a manner, that
        ///  invert_states - ket == bra and vice versa.
        ///  n_bra specifies how many transition graphes there are
        enum state_enum {
              start_state = 0
            , n_bra = S_ORDER
            , n_states = n_bra * 2
            , invert_state = n_states - 1
        };
        ///  \brief contains the information about the bonds
        ///  
        ///  the bonds infront of n_states define what neighbor relations are possible. e.g. down, left, right, up
        ///  would define a square grid. if diag_up and diag_down is added, one gets a triangular grid.
        ///  with just up, down and hori (for horizontal) one gets the hexagonal grid.
        ///  the bonds have to be arranged in such a manner, that
        ///  invert_bond - bond == opposite_bond. e.g. up and down are opposite
        #if GRID_TYPE == 3
        enum bond_enum {
              me = 0
            , start_bond = 1
            , down = 1
            , diag_down
            , right
            , left
            , diag_up
            , up
            , n_bonds
            , hori
            , none
            , middle = left
            , invert_bond = n_bonds - 1 + start_bond
            , tri = 6 + start_bond
            , sqr = 4 + start_bond
            , hex = 3 + start_bond
        };
        #elif GRID_TYPE == 4
        enum bond_enum {
              me = 0
            , start_bond = 1
            , down = 1
            , right
            , left
            , up
            , n_bonds
            , diag_down
            , diag_up
            , hori
            , none
            , middle = left
            , invert_bond = n_bonds - 1 + start_bond
            , tri = 6 + start_bond
            , sqr = 4 + start_bond
            , hex = 3 + start_bond
        };
        #elif GRID_TYPE == 6
        enum bond_enum {
              me = 0
            , start_bond = 1
            , down = 1
            , hori
            , up
            , n_bonds
            , diag_down
            , right
            , left
            , diag_up
            , none
            , middle = hori
            , invert_bond = n_bonds - 1 + start_bond
            , tri = 6 + start_bond
            , sqr = 4 + start_bond
            , hex = 3 + start_bond
        };
        #endif
        ///  \brief contains the information about the spins
        ///  
        ///  the spins have to be arranged in such a manner, that
        ///  invert_spin - spin == opposite_spin
        enum spin_enum {
              start_spin = 0
            , beta = 0
            , alpha
            , n_spins
            , invert_spin = n_spins - 1
        };
        
        enum shift_enum {
              start_shift = 0
            , ket_preswap = 0
            , ket_swap
            , n_shifts
            , no_shift = n_shifts + 1
        };
        
        enum alpha_enum {
              spin_checked = 1
            , all_good = 2
            , bad_spin = 2
            , bad_bond = 4
            , clear = 3
            , not_used = 16
        };
    }

    typedef int spin_type; ///< the spin type, for now just an int
    typedef unsigned loop_type; ///< used for the loop label
    typedef unsigned bond_type; ///< based on bond_enum, but since the enum is not usable as an index, its an size_t
    typedef std::bitset<qmc::n_states> check_type; ///< used for the check variable. a bool would also work, but an size_t8_t needs the same space and provides more option
    typedef unsigned state_type; ///< names the type of the state. again, casting from and to enum all the time would be cumbersome
    typedef unsigned shift_type;
    
    template<typename site_type, int T>
    struct tile_struct {
    };
    //  +---------------------------------------------------+
    //  |            spezialisation for the hex             |
    //  +---------------------------------------------------+
    template<typename site_type>
    struct tile_struct<site_type, qmc::hex>: public std::bitset<6> {
        typedef std::bitset<6> base_type;
    private:
        enum tile_enum_hex {
              bond0 = 0
            , bond1
            , bond2
            , bond3
            , bond4
            , bond5
        };
    public:
        typedef int alpha_type;
        tile_struct(): alpha(0) {
        }
        //------------------- constants -------------------
        static unsigned const tile_per_site = 1;
        static unsigned const n_patterns = 2;
        
        static std::bitset<6> const patterns[2];
        
        void check_bad_bond() {
            for(unsigned i = 0; i < n_patterns; ++i) {
                if((*this) == patterns[i]) {
                    CLEAR_BIT(alpha, qmc::bad_bond)
                    return;
                }
            }
            SET_BIT(alpha, qmc::bad_bond)
        }
        void check_bad_spin() {
            SET_BIT(alpha, qmc::spin_checked)
            
            if(    site->spin[state] != qmc::invert_spin - site->neighbor[qmc::up]->spin[state]
                or site->spin[state] != qmc::invert_spin - site->neighbor[qmc::down]->spin[state]
                or site->spin[state] != qmc::invert_spin - site->neighbor[qmc::down]
                                                               ->neighbor[qmc::hori]
                                                               ->neighbor[qmc::up]
                                                               ->spin[state]
                or site->spin[state] != site->neighbor[qmc::up]->neighbor[qmc::hori]->spin[state]
                or site->spin[state] != site->neighbor[qmc::down]->neighbor[qmc::hori]->spin[state]
            ) {
                SET_BIT(alpha, qmc::bad_spin)
            }
            else
                CLEAR_BIT(alpha, qmc::bad_spin)
        }
        void check_bad_spin_tile(site_type * t, site_type * np, site_type * far) {
            SET_BIT(alpha, qmc::spin_checked)
            if(    t->spin[state] != qmc::invert_spin - np->spin[state]
                or t->spin[state] != far->spin[state])
                SET_BIT(alpha, qmc::bad_spin)
            else
                CLEAR_BIT(alpha, qmc::bad_spin)
        }
        
        bool tile_update() {
            //~ DEBUG_VAR(alpha)
            if(alpha < qmc::all_good) {
                assert(site != NULL);
                if(alpha == 0) {
                    check_bad_spin_tile(site
                                      , site->neighbor[qmc::up + qmc::down - site->bond[state]]
                                      , site->neighbor[site->bond[state]]->neighbor[qmc::hori]
                                       ); //lazy check :-)
                    if(alpha >= qmc::all_good)
                        return false;
                }
                
                //------------------- change bonds -------------------
                site_type * const pos1 = site->neighbor[qmc::up];
                site_type * const pos2 = pos1->neighbor[qmc::hori];
                site_type * const pos3 = pos2->neighbor[qmc::down];
                site_type * const pos4 = pos3->neighbor[qmc::down];
                site_type * const pos5 = pos4->neighbor[qmc::hori];
                
                
                //~ site->bond[state] = base0_     + base1_ -     site->bond[state];
                
                site->bond[state] = qmc::up   + qmc::down - site->bond[state];
                pos1->bond[state] = qmc::hori + qmc::down - pos1->bond[state];
                pos2->bond[state] = qmc::hori + qmc::down - pos2->bond[state];
                pos3->bond[state] = qmc::up   + qmc::down - pos3->bond[state];
                pos4->bond[state] = qmc::up   + qmc::hori - pos4->bond[state];
                pos5->bond[state] = qmc::up   + qmc::hori - pos5->bond[state];
                
                flip();
                //------------------- change neighbor tiles -------------------
                //i
                                                          pos4->tile[state][0].flip(bond0);
                                     pos5->neighbor[qmc::down]->tile[state][0].flip(bond1);
                site->neighbor[qmc::hori]->neighbor[qmc::down]->tile[state][0].flip(bond2);
                  site->neighbor[qmc::hori]->neighbor[qmc::up]->tile[state][0].flip(bond3);
                                       pos1->neighbor[qmc::up]->tile[state][0].flip(bond4);
                                                          pos2->tile[state][0].flip(bond5);
                
                                                          pos4->tile[state][0].check_bad_bond();
                                     pos5->neighbor[qmc::down]->tile[state][0].check_bad_bond();
                site->neighbor[qmc::hori]->neighbor[qmc::down]->tile[state][0].check_bad_bond();
                  site->neighbor[qmc::hori]->neighbor[qmc::up]->tile[state][0].check_bad_bond();
                                       pos1->neighbor[qmc::up]->tile[state][0].check_bad_bond();
                                                          pos2->tile[state][0].check_bad_bond();
                
                return true;
            }
            return false;
        }
        
        void set_info(site_type * const _site, state_type const & _state, unsigned const & _idx) {
            assert(alpha == 0);
            
            state_type bra = (_state < qmc::n_bra ? _state : qmc::invert_state - _state);
            
            if(_site->loop[bra]) {
                alpha = qmc::not_used;
                site = NULL;
                state = 0;
                DEBUG_MSG("nope")
                return;
            }
            DEBUG_MSG("run")
            state = _state;
            site = _site;
            alpha = 0;
            reset();

            
            set(bond0, site->bond[state] == qmc::up);
            set(bond1, site->neighbor[qmc::up  ]->bond[state] == qmc::hori);
            set(bond2, site->neighbor[qmc::up  ]->neighbor[qmc::hori]->bond[state] == qmc::down);
            set(bond3, site->neighbor[qmc::down]->neighbor[qmc::hori]->bond[state] == qmc::up);
            set(bond4, site->neighbor[qmc::down]->bond[state] == qmc::hori);
            set(bond5, site->bond[state] == qmc::down);
            
            check_bad_bond();
            check_bad_spin();
        }
        template<typename Archive>
        void serialize(Archive & ar) {
            ar & *static_cast<base_type*>(this); //upcast
            ar & alpha;
        }
        
        alpha_type alpha;
        site_type * site;
        state_type state;
    };
    template<typename site_type>
    std::bitset<6> const tile_struct<site_type, qmc::hex>::patterns[2] = {(1<<0) + (1<<2) + (1<<4)
                                                                            , (1<<1) + (1<<3) + (1<<5)};
    //  +---------------------------------------------------+
    //  |            spezialisation for the tri             |
    //  +---------------------------------------------------+
    template<typename site_type>
    struct tile_struct<site_type, qmc::tri>: public std::bitset<qmc::n_bonds> {
        typedef int alpha_type;
        typedef std::bitset<qmc::n_bonds> base_type;
        
        tile_struct(): alpha(0) {
        }
        //------------------- constants -------------------
        static unsigned const tile_per_site = 3;
        static unsigned const n_patterns = 2;
        static std::bitset<qmc::n_bonds> const patterns[3][3];
        
        void check_bad_bond() {
            for(unsigned i = 0; i < n_patterns; ++i) {
                if((*this) == patterns[idx][i]) {
                    CLEAR_BIT(alpha, qmc::bad_bond)
                    return;
                }
            }
            SET_BIT(alpha, qmc::bad_bond)
        }
        void check_bad_spin() {
            SET_BIT(alpha, qmc::spin_checked)
            if(    site->spin[state] != qmc::invert_spin - site->neighbor[base1_]->spin[state] 
                or site->spin[state] != qmc::invert_spin - site->neighbor[base0_]->spin[state]
                or site->neighbor[base1_]->spin[state] != qmc::invert_spin - site->neighbor[base1_]->neighbor[base0_]->spin[state]
            ) {
                SET_BIT(alpha, qmc::bad_spin)
            }
            else
                CLEAR_BIT(alpha, qmc::bad_spin)
        }
        void check_bad_spin_tile(site_type * t, site_type * np) {
            SET_BIT(alpha, qmc::spin_checked)
            if(t->spin[state] != qmc::invert_spin - np->spin[state])
                SET_BIT(alpha, qmc::bad_spin)
            else
                CLEAR_BIT(alpha, qmc::bad_spin)
        }
        
        bool tile_update() {
            //~ DEBUG_VAR(alpha)
            if(alpha < qmc::all_good) {
                if(alpha == 0) {
                    check_bad_spin_tile(site, site->neighbor[base0_ + base1_ - site->bond[state]]); //lazy check :-)
                    if(alpha >= qmc::all_good)
                        return false;
                }
                
                //------------------- change bonds -------------------
                site_type * const bas0 = site->neighbor[base0_];
                site_type * const bas1 = site->neighbor[base1_];
                site_type * const diag = site->neighbor[diag_];
                
                
                site->bond[state] = base0_     + base1_     - site->bond[state];
                bas1->bond[state] = base0_     + base1_inv_ - bas1->bond[state];
                bas0->bond[state] = base0_inv_ + base1_     - bas0->bond[state];
                diag->bond[state] = base0_inv_ + base1_inv_ - diag->bond[state];

                flip();
                (*this) &= ~patterns[idx][2]; //mask stuff that has to remain zero
                
                //------------------- change neighbor tiles -------------------
                //i
                
                
                                      bas0->tile[state][idx].flip(base0_);
                                      bas1->tile[state][idx].flip(base1_);
                site->neighbor[base0_inv_]->tile[state][idx].flip(base0_inv_);
                site->neighbor[base1_inv_]->tile[state][idx].flip(base1_inv_);
                
                                      bas0->tile[state][idx].check_bad_bond();
                                      bas1->tile[state][idx].check_bad_bond();
                site->neighbor[base0_inv_]->tile[state][idx].check_bad_bond();
                site->neighbor[base1_inv_]->tile[state][idx].check_bad_bond();
                
                //i+1
                                 bas0->tile[state][ip1].flip(diag_inv_);
                                 diag->tile[state][ip1].flip(diag_);
                bas0->neighbor[diag_]->tile[state][ip1].flip(diag_);
                                 site->tile[state][ip1].flip(diag_inv_);
                
                                 bas0->tile[state][ip1].check_bad_bond();
                                 diag->tile[state][ip1].check_bad_bond();
                bas0->neighbor[diag_]->tile[state][ip1].check_bad_bond();
                                 site->tile[state][ip1].check_bad_bond();
                
                //i+2
                                 bas1->tile[state][ip2].flip(diag_inv_);
                                 diag->tile[state][ip2].flip(diag_);
                bas1->neighbor[diag_]->tile[state][ip2].flip(diag_);
                                 site->tile[state][ip2].flip(diag_inv_);
                
                                 bas1->tile[state][ip2].check_bad_bond();
                                 diag->tile[state][ip2].check_bad_bond();
                bas1->neighbor[diag_]->tile[state][ip2].check_bad_bond();
                                 site->tile[state][ip2].check_bad_bond();
                
                return true;
            }
            return false;
        }
        
        void set_info(site_type * const _site, state_type const & _state, unsigned const & _idx) {
            state = _state;
            site = _site;
            idx = _idx;
            alpha = 0;
            reset();
            if(idx == 0) {
                base0_ = qmc::diag_down;
                base1_ = qmc::up;
                diag_ = qmc::right;
                ip1 = 1;
                ip2 = 2;
            } else if(idx == 1) {
                base0_ = qmc::up;
                base1_ = qmc::left;
                diag_ = qmc::diag_up;
                ip1 = 2;
                ip2 = 0;
            } else if(idx == 2) {
                base0_ = qmc::left;
                base1_ = qmc::diag_down;
                diag_ = qmc::down;
                ip1 = 0;
                ip2 = 1;
            }
            
            base0_inv_ = qmc::invert_bond - base0_;
            base1_inv_ = qmc::invert_bond - base1_;
            diag_inv_ = qmc::invert_bond - diag_;
            
            set(base0_,  site->bond[state] == base1_);
            set(base1_, site->bond[state] == base0_);
            set(base1_inv_,  site->neighbor[base1_]->bond[state] == base0_);
            set(base0_inv_,    site->neighbor[base0_ ]->bond[state] == base1_);
            
            check_bad_bond();
            check_bad_spin();
        }
        
        template<typename Archive>
        void serialize(Archive & ar) {
            ar & *static_cast<base_type*>(this); //upcast
            ar & alpha;
        }
        
        alpha_type alpha;
        site_type * site;
        state_type state;
        unsigned idx;
        unsigned ip1;
        unsigned ip2;
        
        bond_type base0_;
        bond_type base1_;
        bond_type base1_inv_;
        bond_type base0_inv_;
        bond_type diag_;
        bond_type diag_inv_;
    };
    template<typename site_type>
    std::bitset<qmc::n_bonds> const tile_struct<site_type, qmc::tri>::patterns[3][3] = {
                                                                     {(1<<qmc::down) + (1<<qmc::up)
                                                                   ,  (1<<qmc::diag_down) + (1<<qmc::diag_up)
                                                                   ,  (1<<qmc::right) + (1<<qmc::left) + 1}
                                                                   
                                                                   , {(1<<qmc::down) + (1<<qmc::up)
                                                                   ,  (1<<qmc::right) + (1<<qmc::left)
                                                                   ,  (1<<qmc::diag_down) + (1<<qmc::diag_up) + 1}
                                                                    
                                                                   , {(1<<qmc::diag_down) + (1<<qmc::diag_up)
                                                                   ,  (1<<qmc::right) + (1<<qmc::left)
                                                                   ,  (1<<qmc::down) + (1<<qmc::up) + 1}
                                                                   };
    //  +---------------------------------------------------+
    //  |            spezialisation for the sqr             |
    //  +---------------------------------------------------+
    template<typename site_type>
    struct tile_struct<site_type, qmc::sqr>: public std::bitset<qmc::n_bonds> {
        typedef int alpha_type;
        typedef std::bitset<qmc::n_bonds> base_type;
        
        tile_struct(): alpha(0) {
        }
        //------------------- constants -------------------
        static unsigned const tile_per_site = 1;
        static unsigned const n_patterns = 2;
        static std::bitset<qmc::n_bonds> const patterns[2];
        
        void check_bad_bond() {
            for(unsigned i = 0; i < n_patterns; ++i) {
                if((*this) == patterns[i]) {
                    CLEAR_BIT(alpha, qmc::bad_bond)
                    return;
                }
            }
            SET_BIT(alpha, qmc::bad_bond)
        }
        void check_bad_spin() {
            if(    site->spin[state] != qmc::invert_spin - site->neighbor[qmc::right]->spin[state]
                or site->spin[state] != qmc::invert_spin - site->neighbor[qmc::down]->spin[state]
                or site->neighbor[qmc::right]->spin[state] != qmc::invert_spin - site->neighbor[qmc::right]->neighbor[qmc::down]->spin[state]
            ) {
                SET_BIT(alpha, qmc::bad_spin)
            }
            else
                CLEAR_BIT(alpha, qmc::bad_spin)
        }
        void check_bad_spin_tile(site_type * t, site_type * np) {
            if(t->spin[state] != qmc::invert_spin - np->spin[state])
                SET_BIT(alpha, qmc::bad_spin)
            else
                CLEAR_BIT(alpha, qmc::bad_spin)
        }
        
        bool tile_update() {
            if(alpha < qmc::all_good) {
                if(alpha == 0) {
                    check_bad_spin_tile(site, site->neighbor[qmc::down + qmc::right - site->bond[state]]); //lazy check :-)
                    SET_BIT(alpha, qmc::spin_checked)
                    if(alpha >= qmc::all_good)
                        return false;
                }
                
                //------------------- change bonds -------------------
                site_type * const bas0 = site->neighbor[qmc::down];
                site_type * const bas1 = site->neighbor[qmc::right];
                site_type * const diag = bas1->neighbor[qmc::down];
                
                
                site->bond[state] = qmc::down + qmc::right - site->bond[state];
                bas1->bond[state] = qmc::down + qmc::left  - bas1->bond[state];
                bas0->bond[state] = qmc::up   + qmc::right - bas0->bond[state];
                diag->bond[state] = qmc::up   + qmc::left  - diag->bond[state];
                
                flip();
                reset(0);
                
                //------------------- change neighbor tiles -------------------
                for(bond_type b = qmc::start_bond; b < qmc::n_bonds; ++b) {
                    site->neighbor[b]->tile[state][0].flip(b);
                    site->neighbor[b]->tile[state][0].check_bad_bond();
                }
                return true;
            }
            return false;
        }
        
        void set_info(site_type * const _site, state_type const & _state, unsigned const & _idx) {
            //sqr doesn't need _idx since only one tile per site
            state = _state;
            site = _site;
            
            set(qmc::down,  site->bond[state] == qmc::right);
            set(qmc::right, site->bond[state] == qmc::down);
            set(qmc::left,  site->neighbor[qmc::right]->bond[state] == qmc::down);
            set(qmc::up,    site->neighbor[qmc::down ]->bond[state] == qmc::right);
            
            check_bad_bond();
            check_bad_spin();
        }
        
        template<typename Archive>
        void serialize(Archive & ar) {
            ar & *static_cast<base_type*>(this); //upcast
            ar & alpha;
        }
        
        alpha_type alpha;
        site_type * site;
        state_type state;
    };
    template<typename site_type>
    std::bitset<qmc::n_bonds> const tile_struct<site_type, qmc::sqr>::patterns[2]  = {(1<<qmc::down) + (1<<qmc::up)
                                                                                   , (1<<qmc::right) + (1<<qmc::left)};

    struct site_struct;
    
    typedef tile_struct<site_struct, qmc::n_bonds> tile_type;
    ///  \brief the type of the sites in the grid
    struct site_struct {
        
        ///  \brief default constructor
        site_struct(): check(0) {
            for(state_type bra = qmc::start_state; bra != qmc::n_bra; ++bra) {
                loop[bra] = 0;
                spin[bra] = qmc::beta;
                spin[qmc::invert_state - bra] = qmc::beta;
                bond[bra] = qmc::none;
                bond[qmc::invert_state - bra] = qmc::none;
            }
        }
        ///  \brief returns the neightbor of state of (*this)
        ///  
        ///  @param state names the state in which one wants to know the entanglement-parter
        site_struct * partner(bond_type const state) {
            return neighbor[bond[state]];
        }
        site_struct * loop_partner(bond_type & state, bond_type & bra, shift_type const & shift_mode) {
            if(shift_mode == qmc::no_shift) {
                //~ std::cout << "state " << state << "/" << "bond[state] = " << bond[state] << "  n = " << neighbor[bond[state]] << "  states: " << qmc::n_states << std::endl;
                last_dir = bond[state];
                return neighbor[bond[state]];
            }
            
            if(state != bra) {//it's a ket
                state_type new_ket = state - shift_region[shift_mode];
                
                if(new_ket < qmc::n_bra) //lazy boundary for now
                    new_ket += qmc::n_bra;
                
                
                site_struct * partner = neighbor[bond[new_ket]];
                last_dir = bond[new_ket];
                if(partner->shift_region[shift_mode] != shift_region[shift_mode]) {
                    bra += qmc::n_bra - (partner->shift_region[shift_mode] - shift_region[shift_mode]);
                    
                    if(bra >= qmc::n_bra) //lazy boundary for now
                        bra %= qmc::n_bra;
                    
                    //~ std::cout << REDB <<  "ket: jump to " << bra << " loop: " << loop[qmc::invert_state - state] << NONE;
                    state = qmc::invert_state - bra;
                }
                else {
                    //~ std::cout << GREENB <<  "ket: no jump " << bra << " loop: " << loop[qmc::invert_state - state] << NONE;
                }
                return partner;
            }
            else {
                //~ std::cout << "  bra: " << bra << " loop: " << loop[state] << " " << std::endl;
                last_dir = bond[state];
                return neighbor[bond[state]];
            }
        }
        ///  \brief prints the spin of state s12 (12 bc it can be bra or ket) to os
        void print(state_type const & s12 = qmc::start_state, std::ostream & os = std::cout) const {
            
            if(qmc::n_bonds == qmc::tri) {
                os << "(";
                //~ os << tile[s12][0] << "-" << (tile[s12][0].alpha < 2 ? GREENB : RED) << tile[s12][0].alpha << NONE;
                //~ os << "/" << tile[s12][1] << "-" << (tile[s12][1].alpha < 2 ? BLUEB : RED) << tile[s12][1].alpha << NONE;
                //~ os << "/" << tile[s12][2] << "-" << (tile[s12][2].alpha < 2 ? YELLOWB : RED) << tile[s12][2].alpha << NONE;
                os << (tile[s12][0].alpha < 2 ? GREENB : RED) << tile[s12][0].alpha << NONE;
                os << (tile[s12][1].alpha < 2 ? BLUEB : RED) << tile[s12][1].alpha << NONE;
                os << (tile[s12][2].alpha < 2 ? YELLOWB : RED) << tile[s12][2].alpha << NONE;
                os << ")";
            }
            else {
                if(tile[s12][0].alpha != qmc::not_used)
                    os << (tile[s12][0].alpha < 2 ? GREEN : RED) << tile[s12][0].alpha << NONE;
                else
                    os << BLUE << 0 << NONE;
            }
        }
        ///  \brief the fancy print-function used by the grid_class
        std::vector<std::string> const string_print(unsigned const & L, state_type const & s1, unsigned const & what) const {
            std::vector<std::string> res;
            std::stringstream os;
            
            if(qmc::n_bonds == qmc::sqr) {
                os << "  " << print_bond(qmc::up, "|", "", s1, what) << std::left << std::setw(3) << loop[s1]%1000 << std::right << print_bond(qmc::up, "", " ", s1, what);
                res.push_back(os.str()); 
                os.str("");//reset ss
                os << print_bond(qmc::left, "--", "  ", s1, what) << print_spin(s1, what) << print_bond(qmc::right, "---", "   ", s1, what);
                res.push_back(os.str());
                os.str("");//reset ss
                os << "  " << print_bond(qmc::down, "|", " ", s1, what) << " " << "  ";
                res.push_back(os.str());
                return res;
            }
            if(qmc::n_bonds == qmc::tri) {
                os << print_bond(qmc::diag_up, "\\", " ", s1, what) << " " << print_bond(qmc::up, "/", "", s1, what) << std::left << std::setw(3) << loop[s1]%1000 << std::right << print_bond(qmc::up, "", " ", s1, what);
                res.push_back(os.str()); 
                os.str("");//reset ss
                os << print_bond(qmc::left, "--", "  ", s1, what) << print_spin(s1, what) << print_bond(qmc::right, "---", "   ", s1, what);
                res.push_back(os.str());
                os.str("");//reset ss
                os << "  " << print_bond(qmc::down, "/", " ", s1, what) << " " << print_bond(qmc::diag_down, "\\ ", "  ", s1, what);
                res.push_back(os.str());
                return res;
            }
            if(qmc::n_bonds == qmc::hex) {
                
                if(print_alternate%2) {
                    os << "    " << print_bond(qmc::up, "\\", " ", s1, what) << "    ";
                    res.push_back(os.str()); 
                    os.str("");//reset ss
                    os << " " << std::setw(3) << loop[s1]%1000 << " " << print_spin(s1, what) << print_bond(qmc::hori, "---", "   ", s1, what);
                    res.push_back(os.str());
                    os.str("");//reset ss
                    os << "    " << print_bond(qmc::down, "/", " ", s1, what) << "    ";
                    res.push_back(os.str());
                }
                else {
                    os << "   " << print_bond(qmc::up, "/", " ", s1, what) << "     ";
                    res.push_back(os.str()); 
                    os.str("");//reset ss
                    os << print_bond(qmc::hori, "--", "  ", s1, what) << print_spin(s1, what) << " " << std::left << std::setw(3) << loop[s1]%1000 << std::right << "  ";
                    res.push_back(os.str());
                    os.str("");//reset ss
                    os << "   " << print_bond(qmc::down, "\\", " ", s1, what) << "     ";
                    res.push_back(os.str());
                }
                ++print_alternate;
                if(print_alternate % (L + 1) == 0)
                    ++print_alternate;
                return res;
            }
            
            
        }
        
        bool tile_update(state_type const & state, unsigned const & t) {
            return tile[state][t].tile_update();
        }
        template<typename Archive>
        void serialize(Archive & ar) {
            ar & spin;
            ar & bond;
            ar & loop;
            ar & check;
            ar & shift_mode_print;
            ar & print_alternate;
            ar & shift_region;
            ar & tile;
        }
        
        
        
        spin_type spin[qmc::n_states]; ///< looplabel for state
        loop_type loop[qmc::n_bra]; ///< looplabel for each transitiongraph
        bond_type bond[qmc::n_states]; ///< bond-direction for each state
        site_struct * neighbor[qmc::n_bonds]; ///< pointes structure to determine neighbor relations. same for all states
        check_type check; ///< shared by all states
        
        shift_type shift_region[qmc::n_shifts];
        
        static shift_type shift_mode_print;
        static unsigned print_alternate;
        static bond_type last_dir;
        
        tile_type tile[qmc::n_states][tile_type::tile_per_site];
    private:
        ///  \brief plots the bonds in differente colors, depending how the config is
        std::string print_bond(qmc::bond_enum b, std::string go, std::string no, state_type const & s1, unsigned const & what) const {
            std::stringstream res;
            state_type s2 = qmc::invert_state - s1;
            
            std::string ket_color = GREEN;
            std::string trans_color = WHITE;
            if(shift_mode_print != qmc::no_shift) {
                if(shift_region[shift_mode_print] != 0) {
                    ket_color = YELLOWB;
                    trans_color = WHITEB;
                    s2 -= shift_region[shift_mode_print];
                    if(s2 < qmc::n_bra) //lazy boundary for now
                        s2 += qmc::n_bra;
                }
            }
            
            
            if(what == 0) { //print bra and ket
                if(bond[s1] == b and bond[s2] == b)
                    res << trans_color << go << NONE;
                else
                    if(bond[s1] == b)
                        res << MAGENTA << go << NONE;
                    else
                        if(bond[s2] == b)
                            res << ket_color << go << NONE;
                        else
                            res << no;
            }
            else if(what == 1) { //print bra only
                if(bond[s1] == b)
                    res << MAGENTA << go << NONE;
                else
                    res << no;
            }
            else if(what == 2) { //print ket only
                if(bond[s2] == b)
                    res << ket_color << go << NONE;
                else
                    res << no;
            }
                
            return res.str();
        }
        
        std::string print_spin(state_type const & s1, unsigned const & what) const {
            std::stringstream res;
            
            state_type s2 = qmc::invert_state - s1;
            
            if(what == 0) { //print bra and ket
                if(shift_mode_print != qmc::no_shift and shift_region[shift_mode_print] != 0)
                    res << YELLOWB << spin[s1] << NONE;
                else
                    res << BLUEB << spin[s1] << NONE;
                //~ if(spin[s1] != spin[s2])
                    //~ res << YELLOWB << "X" << NONE;
                //~ else
                    //~ res << (spin[s1] == 0 ? BLUEB : REDB) << spin[s1] << NONE;
            }
            else if(what == 1) { //print bra only
                res << (spin[s1] == 0 ? BLUEB : REDB) << spin[s1] << NONE;
            }
            else if(what == 2) { //print ket only
                res << (spin[s2] == 0 ? BLUEB : REDB) << spin[s2] << NONE;
            }
            
            return res.str();
        }
    };
    
    shift_type site_struct::shift_mode_print = qmc::no_shift;
    unsigned site_struct::print_alternate = 1;
    bond_type site_struct::last_dir = 1;
    
    std::ostream & operator<<(std::ostream & os, site_struct const & site) {
        site.print(qmc::start_state, os);
        return os;
    }
}//end namespace perimeter
#endif //__SITE_STRUCT_HEADER
