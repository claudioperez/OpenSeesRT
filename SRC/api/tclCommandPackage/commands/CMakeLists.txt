
target_sources(OpenSeesRT PRIVATE
# Runtime
    "sequential.cpp"
    "runtime.cpp"
    "commands.cpp"
    "nodes.cpp"
    "domain.cpp"
    "recorder.cpp"

# Modeling
    "modeling/model.cpp"
    "modeling/nodes.cpp"
    "modeling/constraint.cpp"
    "modeling/geomTransf.cpp"
    "modeling/rigidLink.cpp"
    "modeling/element.cpp"
    "modeling/region.cpp"
    # sub-element modeling
    "modeling/nDMaterial.cpp"
    "modeling/section.cpp"
    "modeling/uniaxialMaterial.cpp"
    "modeling/uniaxial.cpp"

# Analysis
    "analysis/algorithm.cpp"
    "analysis/integrator.cpp"
    "analysis/analysis.cpp"
    "analysis/numberer.cpp"
    "analysis/ctest.cpp"
    "analysis/solver.cpp"
    "analysis/solver.hpp"

    "modeling/TclUpdateMaterialCommand.cpp"
    "TclUpdateMaterialStageCommand.cpp"

# LOADS & PATTERNS
    "loading/groundMotion.cpp"
    "loading/element_load.cpp"
    "loading/TclPatternCommand.cpp"
    "loading/TclSeriesCommand.cpp"
    "loading/TclSeriesIntegratorCommand.cpp"
    #"domain/pattern/drm/TclPatternCommand.cpp"
)

