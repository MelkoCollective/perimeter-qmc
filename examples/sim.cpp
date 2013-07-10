// Author:  Mario S. Könz <mskoenz@gmx.net>
// Date:    06.05.2013 12:07:06 EDT
// File:    sim.cpp

#include <iostream>
#include <bash_parameter3_msk.hpp>
#include <sim_class.hpp>

using namespace std;
using namespace perimeter;

int main(int argc, char* argv[])
{
    addon::global_seed.set(0);
    
    auto & p = addon::parameter;
    
    p["mult"] = 1;

    p["H"] = 16;
    p["L"] = 16;
    p["shift"] = "shift.txt";
    p["res"] = "results.txt";
    p["timer_dest"] = 1;

    p.read(argc, argv);
    
    p["term"] = p["mult"] * 100000;
    p["sim"] = p["mult"] * 1000000;
    p["H"] = p["L"];
    
    std::string prog_dir = p["prog_dir"];
    
    remove(std::string(prog_dir + std::string(p["res"])).c_str());
    
    p["shift"] = prog_dir + std::string(p["shift"]);
    p["res"] = prog_dir + std::string(p["res"]);
    
    addon::immortal.set_path(p["prog_dir"]);
    
    sim_class sim(p.get());
    
    sim.run();
    
    return 0;
}
