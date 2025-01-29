cmake_minimum_required(VERSION 3.0.0)

project(myLib VERSION 1.0.0)
# set(TARGET_NAME myFlowLib)
# set(CMAKE_CXX_FLAGS "{CMAKE_CXX_FLAGS} -std=gnu++17")
set(CMAKE_EXECUTABLE_SUFFIX ".exe")

# C++ Standard Setting
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

add_executable(field 
main.cpp
)

target_include_directories(field PRIVATE "${PROJECT_SOURCE_DIR}/include")


