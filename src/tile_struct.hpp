// Author:  Mario S. Könz <mskoenz@gmx.net>
// Date:    10.07.2013 16:11:57 EDT
// File:    tile_struct.hpp

#ifndef __TILE_STRUCT_HEADER
#define __TILE_STRUCT_HEADER

#include <enum_typedef.hpp>

#include <bitset>

namespace perimeter_rvb {
    ///  \brief empty general template to represent the tiles
    template<typename site_type, int T>
    struct tile_struct {
    };
    ///  \brief specialisation for the hexagonal grid
    //  +---------------------------------------------------+
    //  |            spezialisation for the hex             |
    //  +---------------------------------------------------+
    template<typename site_type>
    struct tile_struct<site_type, qmc::hex>: public std::bitset<6> { //inherits from bitset
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
        
        ///  \brief the legal bond pattern
        ///  
        ///  0 = no bond, 1 = bond, legal are 010101 or 101010
        static std::bitset<6> const patterns[2];
        
        ///  \brief checks if the bonds allow update
        void check_bad_bond() {
            for(unsigned i = 0; i < n_patterns; ++i) {
                if((*this) == patterns[i]) {
                    CLEAR_BIT(alpha, qmc::bad_bond)
                    return;
                }
            }
            SET_BIT(alpha, qmc::bad_bond)
        }
        ///  \brief checks if the spins allow update
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
        ///  \brief faster spin_check if already known that sites are good (what constrains the spin config)
        void check_bad_spin_tile(site_type * t, site_type * np, site_type * far) {
            SET_BIT(alpha, qmc::spin_checked) //keep in mind that this tile is now checked
            if(    t->spin[state] != qmc::invert_spin - np->spin[state]
                or t->spin[state] != far->spin[state])
                SET_BIT(alpha, qmc::bad_spin)
            else
                CLEAR_BIT(alpha, qmc::bad_spin)
        }
        ///  \brief updates the tile if possible
        ///  
        ///  returns true if success
        bool tile_update() {
            if(alpha < qmc::all_good) {
                assert(site != NULL);
                if(alpha == 0) { //is spin unchecked
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
        ///  \brief a tile needs to know it's site
        ///  
        ///  will tell the tile what state it is, what site it belongs (idx_ only for triangular)
        void set_info(site_type * const _site, state_type const & _state, unsigned const & _idx) {
            assert(alpha == 0);
            
            state_type bra = (_state < qmc::n_bra ? _state : qmc::invert_state - _state);
            
            if(_site->loop[bra]) { //it's crucial that the loops are initialized as they are in order for this to wrok
                alpha = qmc::not_used;
                site = NULL;
                state = 0;
                return;
            }
            
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
        ///  \brief for the checkpoints
        ///  
        ///  this function is used by the serializer to get and set this object
        template<typename Archive>
        void serialize(Archive & ar) {
            ar & *static_cast<base_type*>(this); //upcast to the inherited bitset
            ar & alpha;
        }
        
        alpha_type alpha;
        site_type * site;
        state_type state;
    };
    template<typename site_type>
    std::bitset<6> const tile_struct<site_type, qmc::hex>::patterns[2] = {(1<<0) + (1<<2) + (1<<4)
                                                                            , (1<<1) + (1<<3) + (1<<5)};
    ///  \brief specialisation for the triangular grid
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
        ///  \brief the legal bond pattern
        ///  
        ///  0 = no bond, 1 = bond
        static std::bitset<qmc::n_bonds> const patterns[3][3];
        
        ///  \brief checks if the bonds allow update
        void check_bad_bond() {
            for(unsigned i = 0; i < n_patterns; ++i) {
                if((*this) == patterns[idx][i]) {
                    CLEAR_BIT(alpha, qmc::bad_bond)
                    return;
                }
            }
            SET_BIT(alpha, qmc::bad_bond)
        }
        ///  \brief checks if the spins allow update
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
        ///  \brief faster spin_check if already known that sites are good (what constrains the spin config)
        void check_bad_spin_tile(site_type * t, site_type * np) {
            SET_BIT(alpha, qmc::spin_checked)
            if(t->spin[state] != qmc::invert_spin - np->spin[state])
                SET_BIT(alpha, qmc::bad_spin)
            else
                CLEAR_BIT(alpha, qmc::bad_spin)
        }
        ///  \brief updates the tile if possible
        ///  
        ///  returns true if success
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
        ///  \brief a tile needs to know it's site
        ///  
        ///  will tell the tile what state it is, what site and tile_index it belongs to
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
        ///  \brief for the checkpoints
        ///  
        ///  this function is used by the serializer to get and set this object
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
    ///  \brief specialisation for the square grid
    template<typename site_type>
    struct tile_struct<site_type, qmc::sqr>: public std::bitset<qmc::n_bonds> {
        typedef int alpha_type;
        typedef std::bitset<qmc::n_bonds> base_type;
        
        tile_struct(): alpha(0) {
        }
        //------------------- constants -------------------
        static unsigned const tile_per_site = 1;
        static unsigned const n_patterns = 2;
        ///  \brief the legal bond pattern
        ///  
        ///  0 = no bond, 1 = bond legal are 0101 and 1010
        static std::bitset<qmc::n_bonds> const patterns[2];
        
        ///  \brief checks if the bonds allow update
        void check_bad_bond() {
            for(unsigned i = 0; i < n_patterns; ++i) {
                if((*this) == patterns[i]) {
                    CLEAR_BIT(alpha, qmc::bad_bond)
                    return;
                }
            }
            SET_BIT(alpha, qmc::bad_bond)
        }
        ///  \brief checks if the spins allow update
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
        ///  \brief faster spin_check if already known that sites are good (what constrains the spin config)
        void check_bad_spin_tile(site_type * t, site_type * np) {
            if(t->spin[state] != qmc::invert_spin - np->spin[state])
                SET_BIT(alpha, qmc::bad_spin)
            else
                CLEAR_BIT(alpha, qmc::bad_spin)
        }
        ///  \brief updates the tile if possible
        ///  
        ///  returns true if success
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
        ///  \brief a tile needs to know it's site
        ///  
        ///  will tell the tile what state it is, what site it belongs (idx_ only for triangular)
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
        ///  \brief for the checkpoints
        ///  
        ///  this function is used by the serializer to get and set this object
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

}//end namespace perimeter_rvb
#endif //__TILE_STRUCT_HEADER
