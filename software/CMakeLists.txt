cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(hbt_bes CXX C)
set(GITHOME_MSU $ENV{GITHOME_MSU})
add_custom_target(extern_commonutils COMMAND make WORKING_DIRECTORY ${GITHOME_MSU}/commonutils/software)
add_custom_target(extern_coral COMMAND make WORKING_DIRECTORY ${GITHOME_MSU}/coral/software)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native -O3")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -W -Wall -Wextra -Wmissing-declarations -Wpointer-arith -Wshadow -Wuninitialized -Winit-self -Wundef -Wcast-align -Wformat=2 -Wold-style-cast -Werror=switch -std=c++14 -mfpmath=sse")
message(STATUS "CXX_FLAGS = " ${CMAKE_CXX_FLAGS})

set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
find_package(GSL 1.15 REQUIRED)
find_package(Eigen3 REQUIRED)
include_directories(SYSTEM ${GSL_INCLUDE_DIR} ${EIGEN3_INCLUDE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/include ${GITHOME_MSU}/commonutils/software/include ${GITHOME_MSU}/coral/software/include)
message("Includes: ${GSL_INCLUDE_DIR} ${EIGEN3_INCLUDE_DIR} ${GITHOME_MSU}/commonutils/software/include ${GITHOME_MSU}/coral/software/include" )

file(GLOB hbt_bes_src ${CMAKE_CURRENT_SOURCE_DIR}/src/**.cc)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "lib")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "lib")
add_library(hbt_bes STATIC ${hbt_bes_src})
target_link_libraries(hbt_bes ${GSL_LIBRARY} ${GSL_CBLAS_LIBRARY})
message("Libraries: ${GSL_LIBRARY} ${GSL_CBLAS_LIBRARY}")



