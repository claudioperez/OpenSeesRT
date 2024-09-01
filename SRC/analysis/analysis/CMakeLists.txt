#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================
target_sources(OPS_Analysis
    PRIVATE
      Analysis.cpp 
      DirectIntegrationAnalysis.cpp 
      DomainDecompositionAnalysis.cpp
      DomainUser.cpp 
      EigenAnalysis.cpp
      ResponseSpectrumAnalysis.cpp
      StaticAnalysis.cpp 
      StaticDomainDecompositionAnalysis.cpp 
      SubstructuringAnalysis.cpp    
      TransientAnalysis.cpp
      TransientDomainDecompositionAnalysis.cpp 
      VariableTimeStepDirectIntegrationAnalysis.cpp
    PUBLIC
      Analysis.h 
      DirectIntegrationAnalysis.h 
      DomainDecompositionAnalysis.h
      DomainUser.h 
      EigenAnalysis.h
      ResponseSpectrumAnalysis.h
      StaticAnalysis.h 
      StaticDomainDecompositionAnalysis.h 
      SubstructuringAnalysis.h    
      TransientAnalysis.h
      TransientDomainDecompositionAnalysis.h 
      VariableTimeStepDirectIntegrationAnalysis.h
)

#add_subdirectory(analysis)
