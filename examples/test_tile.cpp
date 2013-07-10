// Author:  Mario S. Könz <mskoenz@gmx.net>
// Date:    03.06.2013 15:55:44 EDT
// File:    test_tile.cpp

#include <iostream>
#include <sim_class.hpp>
#include <bash_parameter3_msk.hpp>
//~ #include <progress_save_msk.hpp>

using namespace std;
using namespace perimeter_rvb;

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
    
    char i, j, s, w;
    char t = 'a';
    
    g.print({0});
    g.print_all({0}, addon::parameter["f"]);    
    
    std::cout << "type: x y state tile cmd" << std::endl;
    std::cout << "possible cmd are: r(two_bond_update), s(invert spin), l(initialize loops), c(update tiles/copy to ket), u(do an update like in a sim), q(quit)" << std::endl;
    
    while(t != 'q') {
        cin >> i >> j >> s >> w >> t;
        i-='0';
        j-='0';
        s-='0';
        w-='0';
        if(t == 'r')
            sim.two_bond_update(i, j, s, w);
        else if(t == 'l')
            g.init_loops();
        else if(t == 'c') {
            g.copy_to_ket();
            g.init_loops();
            g.clear_tile_spin();
        }
        else if(t == 's')
            g(i, j).spin[int(s)] = qmc::invert_spin - g(i, j).spin[int(s)];
        else if(t == 'u')
            sim.update();
        
        g.print({0});
        g.print_all({0}, addon::parameter["f"]);
    }
        
    
    return 0;
}
