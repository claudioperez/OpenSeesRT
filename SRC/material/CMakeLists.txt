#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================
#add_library(OPS_Material OBJECT)

target_sources(OPS_Material
    PRIVATE
      Material.cpp
    PUBLIC
      Material.h
)

target_include_directories(OPS_Material PUBLIC ${CMAKE_CURRENT_LIST_DIR})

target_include_directories(OPS_Material PRIVATE
  ${OPS_SRC_DIR}/material/Solid
  ${OPS_SRC_DIR}/material/nD/cyclicSoil
  ${OPS_SRC_DIR}/material/nD/feap
  ${OPS_SRC_DIR}/material/nD/matCMM
  ${OPS_SRC_DIR}/material/nD/soil
  ${OPS_SRC_DIR}/material/section/fiber
  ${OPS_SRC_DIR}/material/section/yieldSurface
  ${OPS_SRC_DIR}/material/Yield/evolution
  ${OPS_SRC_DIR}/material/Yield/plasticHardeningMaterial
  ${OPS_SRC_DIR}/material/Yield/yieldSurfaceBC
  ${OPS_SRC_DIR}/material/section/integration
  ${OPS_SRC_DIR}/material/section/repres/cell
  ${OPS_SRC_DIR}/material/section/repres/patch
  ${OPS_SRC_DIR}/material/section/repres/reinfBar
  ${OPS_SRC_DIR}/material/section/repres/reinfLayer
  ${OPS_SRC_DIR}/material/section/repres/section
  ${OPS_SRC_DIR}/material/uniaxial/backbone
  ${OPS_SRC_DIR}/material/uniaxial/drain
  ${OPS_SRC_DIR}/material/uniaxial/fedeas
  ${OPS_SRC_DIR}/material/uniaxial/limitState
  ${OPS_SRC_DIR}/material/uniaxial/limitState/limitCurve
  ${OPS_SRC_DIR}/material/uniaxial/PY
  ${OPS_SRC_DIR}/material/uniaxial/snap
  ${OPS_SRC_DIR}/material/uniaxial/stiffness
  ${OPS_SRC_DIR}/material/uniaxial/strength
  ${OPS_SRC_DIR}/material/uniaxial/unloading
  ${OPS_SRC_DIR}/material/nD/stressDensityModel
  ${OPS_SRC_DIR}/material/nD/UANDESmaterials
  ${OPS_SRC_DIR}/material/nD/UWmaterials
)

add_subdirectory(uniaxial)
add_subdirectory(nD)
add_subdirectory(section)
add_subdirectory(Shell)
add_subdirectory(Plane)
add_subdirectory(Yield)

