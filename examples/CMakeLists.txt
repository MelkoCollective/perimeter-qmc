FOREACH (name   sim
                viz
                test_tile
                serialization
                test_serialize
    )
    add_executable(${name} ${name}.cpp)
    #~ add_dependencies(${name} build_nr) #increases the build-number
    add_dependencies(${name} ${SRC}/conf.hpp.in)
    #~ target_link_libraries(${name} ${LIB})
    #~ target_link_libraries(${name} ${LIB} "/usr/lib/libboost_timer.so" "/usr/lib/libboost_system.so.1.50.0")
    add_default_test(${name})
ENDFOREACH(name)
