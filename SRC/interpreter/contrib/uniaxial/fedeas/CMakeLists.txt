#
# FEDEAS uniaxial materials
#
set(fedeas_uniaxial "/home/claudio/packages/FEDEAS/build/src/uniaxial/libfedeas_uniaxial.a")

target_sources(OpenSeesRT PRIVATE
  TclFedeasMaterialCommand.cpp
)
if (EXISTS "${fedeas_uniaxial}")
  target_link_libraries(OpenSeesRT PRIVATE
    OPS_Uniaxial_Fedeas
    ${fedeas_uniaxial}
  )
endif()
