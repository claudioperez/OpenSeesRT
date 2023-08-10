#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================
#                           Select Executable
#
#==============================================================================
set(OPS_FINAL_TARGET "G3" CACHE STRING "OpenSees final target")
#==============================================================================
#                            Basic Switches
#
#==============================================================================

# Optional Extensions
#--------------------------------------
option(OPS_Use_DRM
  "DRM lib"                                                ON )

option(OPS_Use_HDF5
  "HDF5 Dependent Code"                                    OFF)

option(FMK
  "Special FMK Code"                                       OFF)

#==============================================================================
#                           Select Auxiliary Components
#
# Each element in this list owns an associated macro definition that is the
# name of the lib converted to uppercase and prepended with "OPSDEF_"
# (e.g. using OPS_Element_truss defines the macro OPSDEF_ELEMENT_TRUSS)
#==============================================================================

set(OPS_SysOfEqn_List
# OPS_SysOfEqn_UMF
  #OPS_SysOfEqn_ITPACK
)

set(OPS_Extension_List
  OPS_ASDEA
)

# set(OPS_Exclude_List
#   OPS_Element_feap
#   OPS_Material_StressDensity
#   OPS_Recorder_PVD
#   # OPS_Uniaxial_Fedeas
# )

