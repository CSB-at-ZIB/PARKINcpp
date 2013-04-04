ENABLE_LANGUAGE (Fortran)
INCLUDE(FortranCInterface) 

# FIND_PACKAGE(BLAS REQUIRED)
# FIND_PACKAGE(LAPACK REQUIRED)

set(CMAKE_Fortran_FLAGS -fPIC)

add_library(METAN1_A
            metan1.f 
            linalg_metan1.f 
            zibconst.f) 

TARGET_LINK_LIBRARIES(METAN1_A -lgfortran)
# TARGET_LINK_LIBRARIES(METAN1_A ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} -lgfortran)
