project(hbt CXX)
cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
set(GITHOME_MSU $ENV{GITHOME_MSU})

add_custom_target(extern_commonutils COMMAND make WORKING_DIRECTORY ${GITHOME_MSU}/commonutils/software)
add_custom_target(extern_coral COMMAND make WORKING_DIRECTORY ${GITHOME_MSU}/coral/software)
add_custom_target(extern_hbt_bes COMMAND make WORKING_DIRECTORY ${GITHOME_MSU}/ppbar/software)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native -O3 -std=c++17 -fpie -fopenmp")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -W -Wall -Wextra -Wmissing-declarations -Wpointer-arith -Wshadow -Wuninitialized -Winit-self -Wundef -Wcast-align -Wformat=2 -Wold-style-cast -Werror=switch -std=c++11 -mfpmath=sse")
message(STATUS "CXX_FLAGS = " ${CMAKE_CXX_FLAGS})
set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
find_package(GSL 1.15 REQUIRED)
find_package(Eigen3 REQUIRED)
include_directories(SYSTEM ${GSL_INCLUDE_DIR} ${EIGEN3_INCLUDE_DIR} ${GITHOME_MSU}/ppbar/software/include ${GITHOME_MSU}/commonutils/software/include ${GITHOME_MSU}/coral/software/include)          
message("Includes: ${GSL_INCLUDE_DIR} ${EIGEN3_INCLUDE_DIR} ${GITHOME_MSU}/ppbar/software/include ${GITHOME_MSU}/commonutils/software/include ${GITHOME_MSU}/coral/software/include")

link_directories(${GITHOME_MSU}/ppbar/software/lib ${GITHOME_MSU}/commonutils/software/lib ${GITHOME_MSU}/coral/software/lib)
message("Libraries: b3d_bes  ${GSL_LIBRARY} ${GSL_CBLAS_LIBRARY}" )

find_library(${GITHOME_MSU}/commonutils/software/lib commonutils)
find_library(${GITHOME_MSU}/coral/software/lib coral)
find_library(${GITHOME_MSU}/ppbar/software/lib hbt_bes)

#################################################
include_directories(SYSTEM ${GSL_INCLUDE_DIR} ${EIGEN3_INCLUDE_DIR} ${GITHOME_MSU}/ppbar/software/include ${GITHOME_MSU}/commonutils/software/include ${GITHOME_MSU}/coral/software/include)   
link_directories(${GITHOME_MSU}/ppbar/software/lib ${GITHOME_MSU}/commonutils/software/lib ${GITHOME_MSU}/coral/software/lib)

message("Libraries: b3d_bes  ${GSL_LIBRARY} ${GSL_CBLAS_LIBRARY}" )

set(SOURCE_FILES ./psicheck.cc)
SET_SOURCE_FILES_PROPERTIES(
	psicheck.cc PROPERTIES OBJECT_DEPENDS ${GITHOME_MSU}/ppbar/software/lib/libhbt_bes.a OBJECT_DEPENDS ${GITHOME_MSU}/commonutils/software/lib/libcommonutils.a OBJECT_DEPENDS ${GITHOME_MSU}/coral/software/lib/libcoral.a
)

add_executable(psicheck ${SOURCE_FILES})
add_dependencies(psicheck extern_hbt_bes extern_coral extern_commonutils)
target_link_libraries(psicheck hbt_bes coral commonutils ${GSL_LIBRARY} ${GSL_CBLAS_LIBRARY} ${EIGEN3_LIBRARY})
