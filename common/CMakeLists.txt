cmake_minimum_required (VERSION 2.6)

add_library(common
            PARKINLog.cpp)

TARGET_LINK_LIBRARIES(common dlib)