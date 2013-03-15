cmake_minimum_required (VERSION 2.6)

add_library(odelib
            DOP853.cpp
            LIMEXTrajectory.cpp
            LIMEX_A.cpp
            CubicHermiteTrajectory)
            
TARGET_LINK_LIBRARIES(odelib common LIMEX4_3A Ode dlib)