// Author:  Mario S. Könz <mskoenz@gmx.net>
// Date:    06.05.2013 12:07:06 EDT
// File:    sim_playground.cpp

#include <iostream>
#include <sim_class.hpp>
#include <bash_parameter2_msk.hpp>

using namespace std;
using namespace perimeter;

int main(int argc, char* argv[])
{
    std::string test = "";
    //~ std::string test = "../../";
    addon::parameter.set("init0", 0);
    addon::parameter.set("init1", 0);
    addon::parameter.set("p", 0.5);
    addon::parameter.set("f", 1);
    addon::parameter.set("g", 0);
    
    addon::parameter.set("term", 100);
    addon::parameter.set("sim", 1000);
    
    addon::parameter.set("H", 4);
    addon::parameter.set("L", 4);
    
    addon::parameter.read(argc, argv);
    sim_class sim(addon::parameter.get());
    
    sim.grid().print_all({qmc::bra, qmc::bra2}, addon::parameter["-f"]);
    sim.run();
    sim.grid().print_all({qmc::bra, qmc::bra2}, addon::parameter["-f"]);
    sim.grid().print_all({qmc::swap_bra1, qmc::swap_bra2}, addon::parameter["-f"]);
    sim.present_data();
    
    return 0;
}
