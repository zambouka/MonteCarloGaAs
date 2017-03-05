cmake_minimum_required(VERSION 3.6)
project(MCGaAs)

set(CMAKE_CXX_STANDARD 14)

set(SOURCE_FILES main.cpp electron.cpp electron.h mainArea.cpp mainArea.h)
add_executable(MCGaAs ${SOURCE_FILES})