#=====================================================================
# 
#   OpenSees -- Open System For Earthquake Engineering Simulation
#           Pacific Earthquake Engineering Research Center
#
#=====================================================================
target_sources(OpenSeesRT PRIVATE
# Interpreter
    "formats.cpp"
    "interpreter.cpp"
    "pragma.cpp"
    "sequential.cpp"

# Domain
    "commands.cpp"

    "domain/domain.cpp"
    "domain/nodes.cpp"
    "domain/runtime.cpp"
    "domain/printing.cpp"
    "domain/recorder.cpp"
    "domain/TclUpdateMaterialStageCommand.cpp"

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

# Analysis
    "analysis/algorithm.cpp"
    "analysis/integrator.cpp"
    "analysis/analysis.cpp"
    "analysis/numberer.cpp"
    "analysis/ctest.cpp"
    "analysis/solver.cpp"
    "analysis/solver.hpp"

    "modeling/TclUpdateMaterialCommand.cpp"

# LOADS & PATTERNS
    "loading/groundMotion.cpp"
    "loading/element_load.cpp"
    "loading/pattern.cpp"
    "loading/series.cpp"
    "loading/TclSeriesIntegratorCommand.cpp"
    #"domain/pattern/drm/TclPatternCommand.cpp"
)

