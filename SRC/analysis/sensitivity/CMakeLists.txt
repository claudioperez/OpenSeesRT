add_library(OPS_Sensitivity OBJECT)

target_sources(OPS_Sensitivity
    PRIVATE
        NewmarkSensitivityIntegrator.cpp
        SensitivityAlgorithm.cpp
        SensitivityIntegrator.cpp
        StaticSensitivityIntegrator.cpp
    PUBLIC
        NewmarkSensitivityIntegrator.h
        SensitivityAlgorithm.h
        SensitivityIntegrator.h
        StaticSensitivityIntegrator.h
)


target_include_directories(OPS_Sensitivity PUBLIC ${CMAKE_CURRENT_LIST_DIR})


