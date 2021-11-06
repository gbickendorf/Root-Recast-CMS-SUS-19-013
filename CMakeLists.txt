# only for cmake --version >= 3.5.1
cmake_minimum_required(VERSION 3.5.1)
#set(CMAKE_C_COMPILER g++)
#SET (CMAKE_C_COMPILER_WORKS 1)
#SET (CMAKE_CXX_COMPILER_WORKS 1)
# project name
project(main)

# creates the variable EXEC and sets it to hello
set(EXEC main)

# set the C++17 standard
#set(CMAKE_CXX_STANDARD 14)
#set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_C_FLAGS_DEBUG "-g -DDEBUG")
set(CMAKE_CXX_FLAGS_DEBUG "-g -DDEBUG")

#set(CMAKE_CXX_FLAGS "-O3")

# I../includes
find_package(ROOT)

include_directories(include)

add_compile_options(-c -g -Wall )

execute_process(COMMAND "root-config" "--cflags" OUTPUT_VARIABLE ROOT_FLAGS)
separate_arguments(ROOT_FLAGS UNIX_COMMAND "${ROOT_FLAGS}")
add_compile_options(${ROOT_FLAGS})

execute_process(COMMAND "fastjet-config" "--cxxflags" "--libs" "--plugins" OUTPUT_VARIABLE FASTJET_FLAGS)
separate_arguments(FASTJET_FLAGS UNIX_COMMAND "${FASTJET_FLAGS}")
add_compile_options(${FASTJET_FLAGS})

add_compile_options(-llibPhysics -I/home/gerrit/Documents/PHD/Tools/MG5_aMC_v3_2_0/ExRootAnalysis -I/home/gerrit/Documents/PHD/Tools/Delphes-3.5.0)

link_libraries(-L/home/gerrit/root/lib -lGui -lCore -lImt -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lROOTVecOps -lTree -lTreePlayer -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lMultiProc -lROOTDataFrame -pthread -lm -ldl -rdynamic)
link_libraries(-L/home/gerrit/Documents/PHD/Tools/Delphes-3.5.0 -lDelphes -lPhysics)

file(GLOB SOURCES src/*.cpp)

# compiles the files defined by SOURCES to generante the executable defined by EXEC
add_executable(${EXEC} ${SOURCES})

