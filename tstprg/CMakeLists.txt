ENABLE_LANGUAGE (Fortran)

add_executable(tstprg
               main.cpp
               testlinalg.cpp
               testnonlin.cpp
               testparkin_aux.cpp
               testpfizer_simple.cpp
               testsystem.cpp
               testsystem_aux.cpp
               testfoerster_react_c)
               
TARGET_LINK_LIBRARIES(tstprg system nonlin LIMEX4_3A odelib ${CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES})
