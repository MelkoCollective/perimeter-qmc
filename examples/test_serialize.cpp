// Author:  Mario S. Könz <mskoenz@gmx.net>
// Date:    20.06.2013 10:48:49 EDT
// File:    test_serialize.cpp

#include <immortal_msk.hpp>
#include <progress_save_msk.hpp>
#include <bash_parameter3_msk.hpp>
#include <sim_class.hpp>

#include <iostream>
#include <fstream>

int main(int argc, char* argv[]) {
    auto & p = addon::parameter;
    p["init0"] = 0;
    p["init1"] = 0;
    p["f"] = 1;
    p["g"] = 0;

    p["mult"] = 1;

    p["H"] = 4;
    p["L"] = 4;
    p["shift"] = "shift.txt";
    p["spacing"] = 1;
    p["res"] = "results.txt";
    p["timer_dest"] = 1;

    p.read(argc, argv);
    
    p["term"] = p["mult"] * 100000;
    p["sim"] = p["mult"] * 1000000;
    p["H"] = p["L"];
    
    
    addon::checkpoint.set_path(p["prog_dir"]);
    addon::checkpoint.write();
    
    addon::fserial_class backup1("backup_1.bin");
    addon::fserial_class backup2("backup_2.bin");
    
    perimeter_rvb::sim_class sim(p.get());
    perimeter_rvb::sim_class sim2(p.get());
    
    sim.update();
    sim.measure();
    sim.grid().print_all();
    backup1 << sim;
    //~ backup2 << sim;
    
    backup1 >> sim2;
    //~ backup2 >> sim2;
    sim2.grid().print_all();
    
    
    addon::immortal.set_path(p["prog_dir"]);
    
    int data = 1;
    if(addon::immortal.available()) {
        std::cout << "load data" << std::endl;
        addon::immortal >> data;
    }
    char c = 0;
    
    for(unsigned i = addon::immortal.get_index(); i < 100000; ++i) {
        std::cin >> c;
        --data;
        std::cout << i << "  " << data << std::endl;
        //~ if(i%10==9) {
            addon::immortal << data;
            addon::immortal.write_next_index(i + 1);
        //~ }
        if(c == 'q')
            break;
        if(c == 'r') {
            addon::immortal.reset();
            i = 0;
            data = 0;
        }
    }
    
    
    
    
    
    
    return 0;
}
