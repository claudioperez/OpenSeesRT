add_library(UCSD_Array)
target_sources(UCSD_Array
    PRIVATE
        basics.cpp
        BJmatrix.cpp
        BJtensor.cpp
        BJvector.cpp
#       Cosseratstraint.cpp
#       Cosseratstresst.cpp
        nDarray.cpp
        straint.cpp
        stresst.cpp
    PUBLIC
)

#       BJmath_tst.cpp
#       elast_tst.cpp
#       el-pl-test.cpp
#       math_tst.cpp
