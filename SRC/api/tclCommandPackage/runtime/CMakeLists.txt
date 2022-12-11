#=====================================================================
# 
#   OpenSees -- Open System For Earthquake Engineering Simulation
#           Pacific Earthquake Engineering Research Center
#
#=====================================================================

target_sources(OpenSeesRT
    PRIVATE
      BasicAnalysisBuilder.cpp
      BasicModelBuilder.cpp

      # Model Builders
      modelbuilder/TclBuilder.cpp
      modelbuilder/basic/TclBasicBuilder.cpp
      modelbuilder/uniax/TclUniaxialMaterialTester.cpp
      modelbuilder/sect/TclSectionTestBuilder.cpp
      #modelbuilder/plane/TclPlaneStressMaterialTester.cpp


    PUBLIC
        BasicAnalysisBuilder.h
        BasicModelBuilder.h
)

add_subdirectory(Logging)
add_subdirectory(Storage)

target_link_libraries(OpenSeesRT PRIVATE G3_Table)

