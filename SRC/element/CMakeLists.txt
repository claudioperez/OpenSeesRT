#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================
#add_library(OPS_Element OBJECT)
target_include_directories(OPS_Element PUBLIC
  ${OPS_SRC_DIR}/element/absorbentBoundaries
  ${OPS_SRC_DIR}/element/adapter
  ${OPS_SRC_DIR}/element/brick
  ${OPS_SRC_DIR}/element/catenaryCable
  ${OPS_SRC_DIR}/element/CEqElement

  ${OPS_SRC_DIR}/element/community
  ${OPS_SRC_DIR}/element/community/HUelements

# ${OPS_SRC_DIR}/element/beam2d
# ${OPS_SRC_DIR}/element/beam3d
# ${OPS_SRC_DIR}/element/dispBeamColumn
# ${OPS_SRC_DIR}/element/elasticBeamColumn
# ${OPS_SRC_DIR}/element/mixedBeamColumn
  ${OPS_SRC_DIR}/element/Frame/UpdatedLagrange
  ${OPS_SRC_DIR}/element/Frame/GradientInelastic
  ${OPS_SRC_DIR}/element/Frame/ForceInterp
  ${OPS_SRC_DIR}/element/Frame/MixedInterp
  ${OPS_SRC_DIR}/element/Frame/DisplInterp
  ${OPS_SRC_DIR}/element/componentElement
  ${OPS_SRC_DIR}/element/dispBeamColumnInt
# ${OPS_SRC_DIR}/element/nonlinearBeamColumn/matrixutil

  ${OPS_SRC_DIR}/element/bearing/elastomeric
  ${OPS_SRC_DIR}/element/bearing/elastomeric
  ${OPS_SRC_DIR}/element/bearing/friction
  ${OPS_SRC_DIR}/element/bearing/friction/frictionModel

  ${OPS_SRC_DIR}/element/feap
  ${OPS_SRC_DIR}/element/quadrilateral
  ${OPS_SRC_DIR}/element/generic
  ${OPS_SRC_DIR}/element/joint
  ${OPS_SRC_DIR}/element/LHMYS
  ${OPS_SRC_DIR}/element/masonry
  ${OPS_SRC_DIR}/element/mvlem
  ${OPS_SRC_DIR}/element/PML
  ${OPS_SRC_DIR}/element/pyMacro
  ${OPS_SRC_DIR}/element/RockingBC
  ${OPS_SRC_DIR}/element/Shell
  ${OPS_SRC_DIR}/element/surfaceLoad
  ${OPS_SRC_DIR}/element/twoNodeLink

  ${OPS_SRC_DIR}/element/tetrahedron
  ${OPS_SRC_DIR}/element/triangle
  ${OPS_SRC_DIR}/element/truss
)

target_sources(OPS_Element
    PRIVATE
      Element.cpp
      ElementalLoad.cpp
      WrapperElement.cpp
    PUBLIC
      Element.h
      ElementalLoad.h
      WrapperElement.h
)

target_sources(OPS_Utilities
    PRIVATE
      Information.cpp
    PUBLIC
      Information.h
)

target_include_directories(OPS_Element PUBLIC ${CMAKE_CURRENT_LIST_DIR})

add_subdirectory(quadrilateral)
add_subdirectory(twoNodeLink)
add_subdirectory(tetrahedron)
add_subdirectory(triangle)
add_subdirectory(Shell)
add_subdirectory(Brick)


add_subdirectory(zeroLength)

add_subdirectory(truss)
add_subdirectory(catenaryCable)

add_subdirectory(Frame)
# add_subdirectory(elasticBeamColumn)
# add_subdirectory(dispBeamColumn)
# add_subdirectory(forceBeamColumn)
# add_subdirectory(mixedBeamColumn)
# add_subdirectory(gradientInelasticBeamColumn)
# add_subdirectory(updatedLagrangianBeamColumn)

add_subdirectory(dispBeamColumnInt)
add_subdirectory(LHMYS)
add_subdirectory(componentElement)

add_subdirectory(PML)
add_subdirectory(RockingBC)
add_subdirectory(absorbentBoundaries)
add_subdirectory(adapter)

add_subdirectory(bearing/friction)
add_subdirectory(bearing/elastomeric)

add_subdirectory(generic)
add_subdirectory(joint)
add_subdirectory(surfaceLoad)

add_subdirectory(feap)

add_subdirectory(mvlem)
add_subdirectory(pyMacro)
add_subdirectory(masonry)
add_subdirectory(CEqElement)

#add_subdirectory(dmglib)
add_subdirectory(UP-ucsd)
add_subdirectory(community/HUelements)
add_subdirectory(community/XMUelements)
#add_subdirectory(community/UWelements)
# add_subdirectory(PFEMElement)

# TODO: Move out
add_subdirectory(matrixutil)

