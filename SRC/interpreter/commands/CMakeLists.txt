#* ****************************************************************** *#
#*    OpenSees - Open System for Earthquake Engineering Simulation    *#
#*          Pacific Earthquake Engineering Research Center            *#
#* ****************************************************************** *#

target_sources(OpenSeesRT PRIVATE
# Interpreter
    "formats.cpp"
    "interpreter.cpp"
    "pragma.cpp"
    "sequential.cpp"
    "packages.cpp"
    "progress.cpp"

# Domain
    "commands.cpp"

    "domain/domain.cpp"
    "domain/nodes.cpp"
    "domain/runtime.cpp"
    "domain/recorder.cpp"
    "domain/TclUpdateMaterialStageCommand.cpp"
    "domain/TclUpdateMaterialCommand.cpp"

# Modeling
    "modeling/model.cpp"
    "modeling/nodes.cpp"
    "modeling/constraint.cpp"
    "modeling/geomTransf.cpp"
    "modeling/rigidLink.cpp"
    "modeling/element.cpp"
    "modeling/region.cpp"
    "modeling/nDMaterial.cpp"
    "modeling/section.cpp"
    "modeling/uniaxialMaterial.cpp"
    "modeling/uniaxial.cpp"
    "modeling/printing.cpp"

    "invoking/invoke.cpp"
    "invoking/invoke_uniaxial.cpp"
    "invoking/invoke_section.cpp"
    "invoking/invoke_stress.cpp"

# Analysis
    "analysis/algorithm.cpp"
    "analysis/integrator.cpp"
    "analysis/analysis.cpp"
    "analysis/numberer.cpp"
    "analysis/ctest.cpp"
    "analysis/solver.cpp"
    "analysis/solver.hpp"

# LOADS & PATTERNS
    "loading/groundMotion.cpp"
    "loading/element_load.cpp"
    "loading/pattern.cpp"
    "loading/series.cpp"
    "loading/TclSeriesIntegratorCommand.cpp"
    #"domain/pattern/drm/TclPatternCommand.cpp"
)


add_library(OpenSeesRT_Parallel)
target_compile_definitions(OpenSeesRT_Parallel PRIVATE 
  _PARALLEL_INTERPRETERS
  _PARALLEL_PROCESSING
)
target_link_libraries(OpenSeesRT_Parallel PRIVATE OPS_Parallel OPS_Actor)

target_sources(OpenSeesRT_Parallel PRIVATE
  parallel.cpp
)
