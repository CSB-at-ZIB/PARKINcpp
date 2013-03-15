ENABLE_LANGUAGE (Fortran)
INCLUDE(FortranCInterface) 

FIND_PACKAGE(BLAS REQUIRED)
FIND_PACKAGE(LAPACK REQUIRED)

set(CMAKE_Fortran_FLAGS -fPIC)

add_library(LIMEX4_3A
            HERMINE.f 
            LIMD.f 
            LIMDHERM.f 
            LIMDHERM_DD.f 
            sLIMDHERM.f 
            LIMEX4_3A.f 
            LIMEX4_3_Dense.f)

TARGET_LINK_LIBRARIES(LIMEX4_3A ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} -lgfortran)