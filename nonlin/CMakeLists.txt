cmake_minimum_required (VERSION 2.6)

add_library(nonlin
            MultipleShootingGN.cpp
            GaussNewton.cpp
            UserFunc.cpp
            YeOldeParkinCore.cpp)
            
TARGET_LINK_LIBRARIES(nonlin common linalg)
