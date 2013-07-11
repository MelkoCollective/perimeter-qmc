// Author:  Mario S. Könz <mskoenz@gmx.net>
// Date:    10.07.2013 16:07:56 EDT
// File:    enum.hpp

#ifndef __ENUM_HEADER
#define __ENUM_HEADER

#include <conf.hpp>

#include <bitset>
#include <assert.h>

#define DEBUG_VAR(x) std::cout << "\033[1;31m" << "  DEBUG_VAR: " << "\033[0;31m" << #x << " = " << x << "\033[0m" << std::endl;
#define DEBUG_MSG(x) std::cout << "\033[1;31m" << "  DEBUG_MSG: " << "\033[0;31m" << x << "\033[0m" << std::endl;

//it doesn't seem like it sets bits, but the enums are powers of two an thus "bits"
#define SET_BIT(x, y) (x) |= (y);
#define CLEAR_BIT(x, y) (x) &= ~(y);

namespace perimeter_rvb {
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
        ///  \brief contains the different shift commands
        ///  
        ///  add other swap options here if needed. One has to change the shift_files
        ///  and perform a second measurement in order to gain something
        enum shift_enum {
              start_shift = 0
            , ket_preswap = 0
            , ket_swap
            , n_shifts
            , no_shift = n_shifts + 1
        };
        
        ///  \brief used in the tiles
        ///  
        ///  flags set in the alpha variable in the tiles in order to keep track
        ///  what tiles are updateable
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
    typedef std::bitset<qmc::n_states> check_type; ///< used for the check variable. A bool array would also work, but the bitset needs only 1/8 of the space
    typedef unsigned state_type; ///< names the type of the state. again, casting from and to enum all the time would be cumbersome
    typedef state_type shift_type; ///< should be the same as state_type

}//end namespace perimeter_rvb

#endif //__ENUM_HEADER
