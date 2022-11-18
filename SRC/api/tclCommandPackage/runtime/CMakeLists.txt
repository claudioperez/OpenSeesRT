
target_sources(OpenSeesRT
    PRIVATE
        BasicAnalysisBuilder.cpp
        BasicModelBuilder.cpp
    PUBLIC
        BasicAnalysisBuilder.h
        BasicModelBuilder.h
)

add_subdirectory(Logging)
