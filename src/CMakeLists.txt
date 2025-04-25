###
# Copyright (c) 2018  T-Life Research Center, Fudan University, Shanghai, China.
# See the accompanying Manual for the contributors and the way to cite this work.
# Comments and suggestions welcome. Please contact
# Dr. Guanghong Zuo <ghzuo@fudan.edu.cn>
# 
# @Author: Dr. Guanghong Zuo
# @Date: 2017-09-01 16:22:53
# @Last Modified By: Dr. Guanghong Zuo
# @Last Modified Time: 2024-12-25 10:52:32
###

cmake_minimum_required(VERSION 3.14)
project (CVNet)
enable_language(CXX)

## For Compile
##########################################################
SET(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)
SET(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)
LINK_DIRECTORIES(${LIBRARY_OUTPUT_PATH})
SET(CMAKE_VERBOSE_MAKEFILE on)
ADD_DEFINITIONS("-std=c++17 -O3")

## static option
option(STATIC "Build as a static library" OFF)
if(STATIC)
  message(STATUS "Building with STATIC=ON")
  set(CMAKE_FIND_LIBRARY_SUFFIXES .a)
  set(CMAKE_EXE_LINKER_FLAGS "-static -static-libstdc++ -static-libgcc")
endif()

### use CANDA library
if (NOT $ENV{CONDA_PREFIX} STREQUAL "")
  message(STATUS "CURRENT CONDA PREFIX: $ENV{CONDA_PREFIX}")
  include_directories(BEFORE $ENV{CONDA_PREFIX}/include)
  # link_directories(BEFORE $ENV{CONDA_PREFIX}/lib)
endif ()

### Check required LIABRARIES
set(REQUIRED_LIBRARIES ${REQUIRED_LIBRARIES} z)
foreach(LIB ${REQUIRED_LIBRARIES})
  find_library(THE${LIB} ${LIB})
  if(THE${LIB})
    message("-- Found library " ${LIB} ": " ${THE${LIB}})
    LINK_LIBRARIES(${LIB})
  else()
    message(FATAL_ERROR ${LIB} " library not found!")
  endif()
endforeach()

### Check the OPENMP
find_package(OpenMP REQUIRED)
if (OPENMP_FOUND)
  set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  if (APPLE)
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -L/usr/local/lib -lomp")
  endif ()
endif()

## download argparse
include(FetchContent)
FetchContent_Declare(
  argparse
  GIT_REPOSITORY https://github.com/p-ranav/argparse.git
  GIT_TAG        v3.1
)
FetchContent_MakeAvailable(argparse)

## include subdirectories
ADD_SUBDIRECTORY(kit)
ADD_SUBDIRECTORY(cvnet)

