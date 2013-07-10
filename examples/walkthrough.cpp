// Author:  Mario S. Könz <mskoenz@gmx.net>
// Date:    09.07.2013 14:15:41 EDT
// File:    walkthrough.cpp

#include <iostream>
#include <sim_class.hpp>

int main(int argc, char* argv[]) {
    
    perimeter_rvb::grid_class g(6, 6, std::vector<unsigned>(perimeter_rvb::qmc::n_bra, 2));
    
    g.print_all();
    g.print_all({0});
    g.print_all({0}, 7);
    
    g(1,2).spin[0];
    
    
    return 0;
}
