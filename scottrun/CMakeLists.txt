cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(read3dhydro CXX)

add_custom_target(extern_msu_sampler COMMAND make WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/../../best_sampler/software)

#include(${CMAKE_CURRENT_SOURCE_DIR}/../software/CMakeLists.txt)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native -O3")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -W -Wall -Wextra -Wmissing-declarations -Wpointer-arith -Wshadow -Wuninitialized -Winit-self -Wundef -Wcast-align -Wformat=2 -Wold-style-cast -Werror=switch -std=c++11 -mfpmath=sse")
message(STATUS "CXX_FLAGS = " ${CMAKE_CXX_FLAGS})

set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
find_package(GSL 1.15 REQUIRED)
find_package(Eigen3 REQUIRED)
include_directories(SYSTEM ${GSL_INCLUDE_DIR} ${EIGEN3_INCLUDE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/../../best_sampler/software/include ${CMAKE_CURRENT_SOURCE_DIR}/../../best_sampler/software/include)
message("Includes: ${GSL_INCLUDE_DIR} ${EIGEN3_INCLUDE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/../software/include ${CMAKE_CURRENT_SOURCE_DIR}/../../best_sampler/software/include")

link_directories(${CMAKE_CURRENT_SOURCE_DIR}/../software/lib ${CMAKE_CURRENT_SOURCE_DIR}/../../best_sampler/software/lib)
message("Libraries: msu_sampler ${GSL_LIBRARY} ${GSL_CBLAS_LIBRARY}" )

set(SOURCE_FILES ./read3dhydro.cc)
SET_SOURCE_FILES_PROPERTIES(
      read3dhydro.cc PROPERTIES OBJECT_DEPENDS
    ${CMAKE_CURRENT_SOURCE_DIR}/../../best_sampler/software/lib/libmsu_sampler.a
 )
add_executable(read3dhydro ${SOURCE_FILES})

add_dependencies(read3dhydro extern_msu_sampler)
find_library(${CMAKE_CURRENT_SOURCE_DIR}/../software/lib msu_sampler)
target_link_libraries(read3dhydro msu_sampler ${GSL_LIBRARY} ${GSL_CBLAS_LIBRARY} ${EIGEN3_LIBRARY})
