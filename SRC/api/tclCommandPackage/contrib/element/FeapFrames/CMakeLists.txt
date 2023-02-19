#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================

target_sources(OPS_Element
  PRIVATE
    fElement.cpp
    fElmt02.cpp
  PUBLIC
    fElement.h
    fElmt02.h
)

target_sources(OPS_ElementFortran
    PRIVATE
        dummyFeap.f
        fillCommon.f
        getCommon.f
        elmt01.f
        elmt02.f
        elmt03.f
        elmt04.f
        elmt05.f
        elmt06.f
        elmt11.f
)


target_include_directories(OPS_Element PUBLIC ${CMAKE_CURRENT_LIST_DIR})


