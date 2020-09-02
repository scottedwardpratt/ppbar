
cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(hbt_ppbar CXX)
set(GITHOME "/Users/scottpratt/git")

add_custom_target(extern_commonutils COMMAND make WORKING_DIRECTORY ${GITHOME}/commonutils/software)
add_custom_target(extern_coral COMMAND make WORKING_DIRECTORY ${GITHOME}/coral/software)
add_custom_target(extern_hbt_bes COMMAND make WORKING_DIRECTORY ${GITHOME}/ppbar/software)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native -O3")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -W -Wall -Wextra -Wmissing-declarations -Wpointer-arith -Wshadow -Wuninitialized -Winit-self -Wundef -Wcast-align -Wformat=2 -Wold-style-cast -Werror=switch -std=c++11 -mfpmath=sse")
message(STATUS "CXX_FLAGS = " ${CMAKE_CXX_FLAGS})
set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
find_package(GSL 1.15 REQUIRED)
find_package(Eigen3 REQUIRED)
include_directories(SYSTEM ${GSL_INCLUDE_DIR} ${EIGEN3_INCLUDE_DIR} ${GITHOME}/ppbar/software/include ${GITHOME}/commonutils/software/include ${GITHOME}/coral/software/include)          
message("Includes: ${GSL_INCLUDE_DIR} ${EIGEN3_INCLUDE_DIR} ${GITHOME}/ppbar/software/include ${GITHOME}/commonutils/software/include ${GITHOME}/coral/software/include")

link_directories(${GITHOME}/ppbar/software/lib ${GITHOME}/commonutils/software/lib ${GITHOME}/coral/software/lib)
message("Libraries: b3d_bes                                     ${GSL_LIBRARY} ${GSL_CBLAS_LIBRARY}" )

find_library(${GITHOME}/commonutils/software/lib commonutils)
find_library(${GITHOME}/coral/software/lib coral)
find_library(../../software/lib hbt_bes)

set(SOURCE_FILES ./hbt_ppbar.cc)
SET_SOURCE_FILES_PROPERTIES(
	hbt_ppbar.cc PROPERTIES OBJECT_DEPENDS
	${GITHOME}/ppbar/software/lib/libhbt_bes.a
)

add_executable(hbt_ppbar ${SOURCE_FILES})
add_dependencies(hbt_ppbar extern_hbt_bes extern_coral extern_commonutils)
target_link_libraries(hbt_ppbar hbt_bes commonutils coral ${GSL_LIBRARY} ${GSL_CBLAS_LIBRARY} ${EIGEN3_LIBRARY})

