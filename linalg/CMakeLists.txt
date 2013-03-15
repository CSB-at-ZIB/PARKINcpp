cmake_minimum_required (VERSION 2.6)

add_library(linalg
            Matrix.cpp
            QRCholesky.cpp
            QRDecomp.cpp
            QRMoorePenrose.cpp
            QRPetersWilkinson.cpp
            QRconDecomp.cpp
            Vector.cpp)
            
TARGET_LINK_LIBRARIES(linalg dlib common)