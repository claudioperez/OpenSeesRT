#=====================================================================
# 
#   OpenSees -- Open System For Earthquake Engineering Simulation
#           Pacific Earthquake Engineering Research Center
#
#=====================================================================

target_sources(OPS_Runtime
    PRIVATE
      TclPackageClassBroker.cpp
      # Model Builders
      BasicAnalysisBuilder.cpp
      BasicModelBuilder.cpp
      modelbuilder/TclBuilder.cpp
      modelbuilder/basic/TclBasicBuilder.cpp

    PUBLIC
      BasicAnalysisBuilder.h
      BasicModelBuilder.h
      TclPackageClassBroker.h
)

add_subdirectory(Logging)
add_subdirectory(Storage)

target_link_libraries(OPS_Runtime PRIVATE G3_Table)

