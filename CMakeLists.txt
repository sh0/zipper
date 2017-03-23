#
# ZipPack - Mesh Zippering
#

# CMake requirements
cmake_minimum_required(VERSION 3.1)

# Paths
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake" "${CMAKE_SOURCE_DIR}/cmake/Modules")
set(SOURCE_DIR "${CMAKE_SOURCE_DIR}/source")

# Project
project(zipper)

# Options
option(BUILD_STATIC "Build static library" OFF)

# Tools
include(tools-compilers)

# Put libraries and executables in same directory
set(EXECUTABLE_OUTPUT_PATH "${PROJECT_BINARY_DIR}/bin")
set(LIBRARY_OUTPUT_PATH "${PROJECT_BINARY_DIR}/bin")

# Subdirs
add_subdirectory(meschach)
add_subdirectory(source)
