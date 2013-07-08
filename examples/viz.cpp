// Author:  Mario S. Könz <mskoenz@gmx.net>
// Date:    13.06.2013 12:06:28 EDT
// File:    vis.cpp

#include <iostream>
#include <progress_save_msk.hpp>
#include <bash_parameter3_msk.hpp>
#include <sim_class.hpp>
#include <stdio.h>

using namespace std;
using namespace perimeter;

int main(int argc, char* argv[])
{
    //~ addon::global_seed.set(0);
    auto & p = addon::parameter;
    p["init0"] = 2;
    p["init1"] = 2;
    p["f"] = 1;
    p["g"] = 0;

    p["mult"] = 1;

    p["H"] = 6;
    p["L"] = 6;
    p["shift"] = "shift.txt";
    p["spaceing"] = 1;
    p["res"] = "results.txt";
    p["timer_dest"] = 0;

    p.read(argc, argv);
    
    p["term"] = p["mult"] * 100000;
    p["sim"] = p["mult"] * 1000000;
    
    std::string prog_dir = p["prog_dir"];
    p["shift"] = prog_dir + std::string(p["shift"]);
    p["res"] = prog_dir + std::string(p["res"]);
    
    sim_class sim(p.get());
    
    remove((prog_dir + "/../../../SimuViz/Example/viz.txt").c_str());
    
    for(size_t i = 0; i < 100; ++i) {
        sim.update();
    }
        
    return 0;
}
