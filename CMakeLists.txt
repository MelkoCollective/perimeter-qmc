cmake_minimum_required (VERSION 2.6)
project(MasterThesis)

#=================== compiler & flags ===================
SET(CMAKE_CXX_COMPILER "g++")
#~ SET(CMAKE_CXX_COMPILER "g++-4.7") #sharcnet?
#~ SET(CMAKE_CXX_COMPILER "cc") #sharcnet
SET(CMAKE_CXX_FLAGS "-Wall -std=c++11 -O3 -pedantic")
#~ SET(CMAKE_EXE_LINKER_FLAGS "-pg") #for profiling

#=================== set paths ===================
set(SRC "${PROJECT_SOURCE_DIR}/src")
set(DOCU "${PROJECT_SOURCE_DIR}/doc")
set(TEST "${PROJECT_SOURCE_DIR}/test")
set(EXAMPLE "${PROJECT_SOURCE_DIR}/test")


set(LIB demo_lib)

set(MAJOR 1)
set(MINOR 1)
#=================== set options ===================
OPTION(BUILD_TEST "build all the tests" ON)
OPTION(BUILD_EXAMPLE "build all the examples" ON)
SET(USE_GRID 3 CACHE STRING "choose the grid type (3=tri, 4=sqr, 6=hex)")
SET(USE_S 2 CACHE STRING "choose the order of the Renyi entropy")


#=================== custom stuff ===================
add_custom_target(
            doc ALL
            COMMAND doxygen "${SRC}/Doxyfile"
            )
            
add_custom_target(
            opendoc ALL
            COMMAND firefox "${DOCU}/html/annotated.html" &
            )

#~ add_custom_target(
            #~ stats ALL
            #~ COMMAND gitstats "../ ../doc/" &
            #~ )

macro(add_default_test)
    set(name ${ARGV0})
    add_test(${name} ${name})
endmacro(add_default_test)

#=================== write config files ===================
configure_file(
"${SRC}/conf.hpp.in"
"${SRC}/conf.hpp"
)
#=================== include directories ===================
include_directories("${SRC}")
include_directories("${SRC}/addons")
include_directories("${SRC}/version")

#=================== added subdirs ===================
add_subdirectory("src")

if(BUILD_TEST)
    #~ enable_testing()
    add_subdirectory("test")
endif(BUILD_TEST)

if(BUILD_EXAMPLE)
    add_subdirectory("examples")
endif(BUILD_EXAMPLE)


