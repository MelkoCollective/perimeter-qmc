FOREACH (name   test_1
    )
    add_executable(${name} ${name}.cpp)
    add_dependencies(${name} ${SRC}/conf.hpp.in)
    add_default_test(${name})
ENDFOREACH(name)
