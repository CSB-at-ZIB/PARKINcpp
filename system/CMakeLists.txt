cmake_minimum_required (VERSION 2.6)

add_library(system
            BioMedication.cpp
            BioPAR.cpp
            BioProcessor.cpp
            BioRHS.cpp
            BioSystem.cpp
            BioSystemODE.cpp
            BioSystemVAR.cpp
            ExprNode.cpp
            Expression.cpp)
            
TARGET_LINK_LIBRARIES(system common linalg nonlin dlib odelib)
