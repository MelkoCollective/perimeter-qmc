// Author:  Mario S. Könz <mskoenz@gmx.net>
// Date:    03.06.2013 15:55:44 EDT
// File:    test_tile.cpp

#include <iostream>
#include <sim_class.hpp>
#include <bash_parameter3_msk.hpp>
//~ #include <progress_save_msk.hpp>

using namespace std;
using namespace perimeter;

int main(int argc, char* argv[])
{
    addon::parameter.set("init0", 0);
    addon::parameter.set("init1", 0);
    addon::parameter.set("f", 7);
    addon::parameter.set("g", 0);

    addon::parameter.set("mult", 1);

    addon::parameter.set("H", 4);
    addon::parameter.set("L", 4);
    addon::parameter.set("shift", "shift.txt");
    addon::parameter.set("res_file", "results.txt");
    
    
    addon::parameter.read(argc, argv);
    
    addon::parameter.set("term", addon::parameter["mult"] * 100000);
    addon::parameter.set("sim", addon::parameter["mult"] * 1000000);
    addon::parameter["shift"] = std::string(addon::parameter["prog_dir"]) + std::string(addon::parameter["shift"]);
    std::cout << addon::parameter["shift"] << std::endl;
    
    
    sim_class sim(addon::parameter.get());
    grid_class & g(sim.grid());
    //~ g(1, 0).spin[0] = 0;
    //~ g(1, 1).spin[0] = 1;
    
    //~ std::cout << g(1, 1).tile_update(0, 1) << std::endl;
    //~ std::cout << g(1, 3).tile_update(0, 1) << std::endl;
    //~ std::cout << g(1, 0).tile_update(0, 1) << std::endl;
    //~ std::cout << g(1, 2).tile_update(0, 1) << std::endl;
    
    char i, j, s, t, w;
    
    while(i!=10) {
        cin >> i >> j >> s >> w >> t;
        i-='0';
        j-='0';
        s-='0';
        w-='0';
        t-='0';
        if(t == 'r'-'0')
            sim.two_bond_update(i, j, s, w);
        else if(t == 'l'-'0')
            g.init_loops();
        else if(t == 'c'-'0') {
            g.copy_to_ket();
            g.init_loops();
            g.clear_tile_spin();
        }
        else if(t == 's'-'0')
            g(i, j).spin[s] = qmc::invert_spin - g(i, j).spin[s];
        else if(t == 'u'-'0')
            sim.update();
        //~ else
            //~ g(i, j).tile_update(0, 0);
            //~ g(i, j).tile_update(0, t);
        
        g.print({0});
        g.print_all({0}, addon::parameter["f"]);
    }
    
    //~ g(2, 3).spin[0] = 0;
    //~ g(3, 3).spin[0] = 1;
    
    //~ g.init_loops();
    
    g.print({0});
    g.print_all({0}, addon::parameter["f"]);    
    
    return 0;
}
