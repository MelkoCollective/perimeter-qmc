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

#include <tile_struct.hpp>

//perimeter is documented in grid_class.hpp
namespace perimeter_rvb {
    
    ///  \brief the representation of the sites
    ///  
    ///  this is actually a super-site that hold sites for all states with the same x and y
    ///  the final 2D array (x, y) of super-site better than to have a 3D array (x, y, state) of actual sites because of the shift
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
        ///  
        ///  this function ignores any shift and just returns the partner in the same state
        site_struct * partner(bond_type const state) {
            return neighbor[bond[state]];
        }
        ///  \brief returns the neightbor of state of (*this) with shift
        ///  
        ///  @param state names the state in which one wants to know the entanglement-parter. May be changed by this function
        ///  @param bra is the corresponding bra to state. If state is already a "bra" then bra == state. May be changed by this function
        ///  @param shift_mode says what shift_mode (preswap/swap) is active
        ///  
        ///  this function returns the partner taking the shift into account
        site_struct * loop_partner(bond_type & state, bond_type & bra, shift_type const & shift_mode) {
            if(shift_mode == qmc::no_shift) {
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
                    
                    state = qmc::invert_state - bra;
                }
                return partner;
            }
            else {
                last_dir = bond[state];
                return neighbor[bond[state]];
            }
        }
        ///  \brief prints the alpha-variable of the tiles in a nice way
        void print(state_type const & s12 = qmc::start_state, std::ostream & os = std::cout) const {
            
            if(qmc::n_bonds == qmc::tri) {
                os << "(";
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
        ///  \brief just forwards to the tile_update function
        ///  
        ///  state and t_nr just specify what tile should get updated
        bool tile_update(state_type const & state, unsigned const & t_nr) {
            return tile[state][t_nr].tile_update();
        }
        ///  \brief for the checkpoints
        ///  
        ///  this function is used by the serializer to get and set this object
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
        check_type check; ///< shared by all states. A bitset, can be used like a vector of bools
        
        shift_type shift_region[qmc::n_shifts]; ///< says by how much the state has to be permuted for the various shift_modes
        tile_type tile[qmc::n_states][tile_type::tile_per_site]; ///< the tiles that are managed by this site
        
        static shift_type shift_mode_print; ///< for nicer printing
        static unsigned print_alternate; ///< for nicer printing
        static bond_type last_dir; ///< for tracking the sign. Shows the direction of the last loop_partner return
        
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
        ///  \brief nicer printing
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
}//end namespace perimeter_rvb
#endif //__SITE_STRUCT_HEADER
