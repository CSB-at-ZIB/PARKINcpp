cmake_minimum_required (VERSION 2.6)

add_library(odelib
            DOP853.cpp
            LIMEXTrajectory.cpp
            LIMEX_A.cpp
            METAN_A.cpp
            LinearTrajectory.cpp
            CubicHermiteTrajectory.cpp)
            
TARGET_LINK_LIBRARIES(odelib common METAN1_A LIMEX4_3A Ode dlib)
