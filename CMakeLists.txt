cmake_minimum_required (VERSION 2.6)

project (PARKINcpp)

if(UNIX)
    set(CMAKE_C_FLAGS   "${CMAKE_C_FLAGS}   -fPIC")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
endif(UNIX)

INCLUDE_DIRECTORIES ("${PROJECT_SOURCE_DIR}")

add_subdirectory (addpkg)
add_subdirectory (common)
add_subdirectory (linalg)
add_subdirectory (nonlin)
add_subdirectory (odelib)
add_subdirectory (system)
add_subdirectory (swig)
add_subdirectory (tstprg)
