
# target_link_libraries(OpenSeesCommunity PRIVATE OPS_Material)

target_include_directories(OpenSeesRT PUBLIC
  ${OPS_SRC_DIR}/material/uniaxial/stiffness
  ${OPS_SRC_DIR}/material/uniaxial/strength
  ${OPS_SRC_DIR}/material/uniaxial/unloading
  ${OPS_SRC_DIR}/material/uniaxial/backbone
)

target_sources(OpenSeesRT PRIVATE
    "backbone.cpp"
    "StiffnessDegradationCommand.cpp"
    "StrengthDegradationCommand.cpp"
    "UnloadingRuleCommand.cpp"
)
