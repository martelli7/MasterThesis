# $Id: CMakeLists.txt 1 2024-12-02 martelli $

#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 3.10 FATAL_ERROR)
project(RootAnalysis)
set(CMAKE_BUILD_TYPE Debug)
#----------------------------------------------------------------------------
# Find ROOT (required package)
#
find_package(ROOT REQUIRED COMPONENTS RIO Net Hist Tree MathCore GenVector)

#----------------------------------------------------------------------------

#include_directories(${ROOT_INCLUDE_DIRS})

include_directories(${PROJECT_SOURCE_DIR}/include)


#----------------------------------------------------------------------------
# Locate sources and headers for this project
# NB: headers are included so they will show up in IDEs
#
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(RootAnalysis RootAnalysis.cc ${sources})
target_link_libraries(RootAnalysis ${ROOT_LIBRARIES})

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS RootAnalysis DESTINATION bin)


