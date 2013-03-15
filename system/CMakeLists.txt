cmake_minimum_required (VERSION 2.6)

add_library(system
            BioPAR.cpp
            BioProcessor.cpp
            BioRHS.cpp
            BioSystem.cpp
            ExprNode.cpp
            Expression.cpp)
            
TARGET_LINK_LIBRARIES(system common linalg nonlin dlib odelib)