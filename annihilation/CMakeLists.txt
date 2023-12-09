project(ppbartest CXX)
cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
set(GITHOME_MSU $ENV{GITHOME_MSU})

add_custom_target(extern_commonutils COMMAND make WORKING_DIRECTORY ${GITHOME_MSU}/commonutils/software)
add_custom_target(extern_coral COMMAND make WORKING_DIRECTORY ${GITHOME_MSU}/coral/software)

set(CMAKE_CXX_FLAGS "-W -Wall -march=native -std=c++17 -O2")
message(STATUS "CXX_FLAGS = " ${CMAKE_CXX_FLAGS})
set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
find_package(GSL REQUIRED)
find_package(Eigen3 REQUIRED)

include_directories(SYSTEM ${GSL_INCLUDE_DIR} ${EIGEN3_INCLUDE_DIR} ${GITHOME_MSU}/ppbar/software/include ${GITHOME_MSU}/commonutils/software/include ${GITHOME_MSU}/coral/software/include)

link_directories(${GITHOME_MSU}/ppbar/software/lib ${GITHOME_MSU}/commonutils/software/lib ${GITHOME_MSU}/coral/software/lib)

find_library(${GITHOME_MSU}/commonutils/software/lib commonutils)
find_library(${GITHOME_MSU}/coral/software/lib coral)

#################################################

link_directories(${GITHOME_MSU}/ppbar/software/lib ${GITHOME_MSU}/commonutils/software/lib ${GITHOME_MSU}/coral/software/lib)

#################################################

set(SOURCE_FILES ./ppbartest.cc)
SET_SOURCE_FILES_PROPERTIES(ppbartest.cc
	PROPERTIES OBJECT_DEPENDS ${GITHOME_MSU}/coral/software/lib/libcoral.a
)
SET_SOURCE_FILES_PROPERTIES(ppbartest.cc
	PROPERTIES OBJECT_DEPENDS  ${GITHOME_MSU}/commonutils/software/lib/libcommonutils.a
)

add_executable(ppbartest ${SOURCE_FILES})
add_dependencies(ppbartest extern_coral extern_commonutils)
target_link_libraries(ppbartest ${GITHOME_MSU}/coral/software/lib/libcoral.a ${GITHOME_MSU}/commonutils/software/lib/libcommonutils.a  ${GSL_LIBRARY} ${GSL_CBLAS_LIBRARY} ${EIGEN3_LIBRARY})

#################################################

set(SOURCE_FILES ./cwcomplex.cc)
SET_SOURCE_FILES_PROPERTIES(cwcomplex.cc
	PROPERTIES OBJECT_DEPENDS ${GITHOME_MSU}/coral/software/lib/libcoral.a
)
SET_SOURCE_FILES_PROPERTIES(cwcomplex.cc
	PROPERTIES OBJECT_DEPENDS  ${GITHOME_MSU}/commonutils/software/lib/libcommonutils.a
)
project(cwcomplex CXX)
add_executable(cwcomplex ${SOURCE_FILES})
add_dependencies(cwcomplex extern_coral extern_commonutils)
target_link_libraries(cwcomplex ${GITHOME_MSU}/coral/software/lib/libcoral.a ${GITHOME_MSU}/commonutils/software/lib/libcommonutils.a  ${GSL_LIBRARY} ${GSL_CBLAS_LIBRARY} ${EIGEN3_LIBRARY})

#################################################

set(SOURCE_FILES ./wftest.cc)
SET_SOURCE_FILES_PROPERTIES(wftest.cc
	PROPERTIES OBJECT_DEPENDS ${GITHOME_MSU}/coral/software/lib/libcoral.a
)
SET_SOURCE_FILES_PROPERTIES(wftest.cc
	PROPERTIES OBJECT_DEPENDS  ${GITHOME_MSU}/commonutils/software/lib/libcommonutils.a
)
project(wftest CXX)
add_executable(wftest ${SOURCE_FILES})
add_dependencies(wftest extern_coral extern_commonutils)
target_link_libraries(wftest ${GITHOME_MSU}/coral/software/lib/libcoral.a ${GITHOME_MSU}/commonutils/software/lib/libcommonutils.a  ${GSL_LIBRARY} ${GSL_CBLAS_LIBRARY} ${EIGEN3_LIBRARY})

#################################################

set(SOURCE_FILES ./coulfit.cc)
SET_SOURCE_FILES_PROPERTIES(coulfit.cc
	PROPERTIES OBJECT_DEPENDS ${GITHOME_MSU}/coral/software/lib/libcoral.a
)
SET_SOURCE_FILES_PROPERTIES(coulfit.cc
	PROPERTIES OBJECT_DEPENDS  ${GITHOME_MSU}/commonutils/software/lib/libcommonutils.a
)
project(coulfit CXX)
add_executable(coulfit ${SOURCE_FILES})
add_dependencies(coulfit extern_coral extern_commonutils)
target_link_libraries(coulfit ${GITHOME_MSU}/coral/software/lib/libcoral.a ${GITHOME_MSU}/commonutils/software/lib/libcommonutils.a  ${GSL_LIBRARY} ${GSL_CBLAS_LIBRARY} ${EIGEN3_LIBRARY})