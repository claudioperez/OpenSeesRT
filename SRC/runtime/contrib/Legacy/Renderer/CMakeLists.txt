#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================
add_library(OPS_Renderer           OBJECT EXCLUDE_FROM_ALL)
add_library(OPS_Graphics           OBJECT EXCLUDE_FROM_ALL)
add_library(OPS_Graphics_Default   OBJECT EXCLUDE_FROM_ALL)
add_library(OPS_Graphics_GL        OBJECT EXCLUDE_FROM_ALL)

target_sources(OPS_Renderer
  PRIVATE
    PlainMap.cpp
    Renderer.cpp
  PUBLIC
    PlainMap.h
    Renderer.h
    #ColorMap.h
)

target_sources(OPS_Graphics
    PRIVATE
      ${OPS_SRC_DIR}/recorder/AlgorithmIncrements.cpp
      ${OPS_SRC_DIR}/recorder/FilePlotter.cpp
    PUBLIC
      ${OPS_SRC_DIR}/recorder/AlgorithmIncrements.h
      ${OPS_SRC_DIR}/recorder/FilePlotter.h
)


target_sources(OPS_Graphics
  PRIVATE
    DofColorMap.cpp
  PUBLIC
    DofColorMap.h
)

target_sources(OPS_Graphics_Default
  PRIVATE
    Clipping.cpp
    db.cpp
    Device.cpp
    gMatrix.cpp
    Projection.cpp
    Scan.h
    View.h
    Viewport.cpp
    VrmlViewer.cpp
    WindowDevice.cpp
    WindowRenderer.cpp
    X11Device.cpp
    X11Renderer.cpp
  PUBLIC
    Clipping.h
    Device.h
    Projection.h
    Scan.cpp
    View.cpp
    Viewport.h
    VrmlViewer.h
    WindowDevice.h
    WindowRenderer.h
    X11Device.h
    X11Renderer.h
)

target_sources(OPS_Graphics_GL
  PRIVATE
    AGL_Device.cpp
    OpenGlDevice.cpp
    OpenGlRenderer.cpp
  PUBLIC
    AGL_Device.h
    OpenGlDevice.h
    OpenGLRenderer.h
)

target_include_directories(OPS_Renderer PUBLIC ${CMAKE_CURRENT_LIST_DIR})

#----------------------------
if (TRUE)
  add_compile_definitions(_NOGRAPHICS)

elseif (${OPS_Use_Graphics_Option} STREQUAL "Base")
  target_link_libraries(${OPS_FINAL_TARGET} OPS_Graphics_Default OPS_Graphics)

elseif (${OPS_Use_Graphics_Option} STREQUAL "OpenGL")
  message("OPS >>> Including OpenGL graphics option")
  set_source_files_properties(
    "${OPS_SRC_DIR}/recorder/AlgorithmIncrements.cpp"
    "${OPS_SRC_DIR}/recorder/FilePlotter.cpp"
    "${OPS_SRC_DIR}/renderer/main.cpp"
    "${OPS_SRC_DIR}/renderer/OpenGlDevice.cpp"
    "${OPS_SRC_DIR}/renderer/OpenGlRenderer.cpp"
    "${OPS_SRC_DIR}/tcl/TclFeViewer.cpp"
    "${OPS_SRC_DIR}/tcl/TclVideoPlayer.cpp"
    PROPERTIES COMPILE_DEFINITIONS _AGL
  )
  target_link_libraries(${OPS_FINAL_TARGET} OPS_Graphics_GL OPS_Graphics)
endif()

