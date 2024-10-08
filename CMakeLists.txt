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
find_package(OpenMP REQUIRED)

find_package(GSL REQUIRED)

add_compile_options(-fopenmp -I/home/gerrit/Documents/PHD/Tools/MG5_aMC_v3_3_2/ExRootAnalysis -I/home/gerrit/Documents/PHD/Tools/Delphes-3.5.0 -I/home/gerrit/Documents/PHD/Tools/boost_1_79_0)

link_libraries(-L/home/gerrit/root/lib -lFoam -lRooFitCore -lRooStats -lRooFit -lHtml -lMinuit -lGui -lCore -lImt -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lROOTVecOps -lTree -lTreePlayer -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lMultiProc -lROOTDataFrame -pthread -lm -ldl -rdynamic)
link_libraries(-L/home/gerrit/Documents/PHD/Tools/Delphes-3.5.0 -lDelphes -lPhysics "${OpenMP_CXX_FLAGS}")
link_libraries(-L/home/gerrit/Documents/PHD/Tools/MG5_aMC_v3_3_2/ExRootAnalysis -lExRootAnalysis)
link_libraries(GSL::gsl GSL::gslcblas)

file(GLOB SOURCES src/*.cpp)

# compiles the files defined by SOURCES to generante the executable defined by EXEC
add_executable(${EXEC} ${SOURCES})

