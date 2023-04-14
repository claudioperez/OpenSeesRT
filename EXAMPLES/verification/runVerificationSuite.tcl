# script to run all the verification scripts
# PASS/FAILURE results in file README.md when run

# open results file such that it is cleared out of any data
set results [open README.md w]
puts $results "| Status | Notes |\n|--------|------------------------------|"
close $results

# cd Basic
# source sdofTransient.tcl
# source SmallEigen.tcl
# source NewmarkIntegrator.tcl
# cd ..

source PlanarTruss.tcl
source PlanarTruss.Extra.tcl
source PortalFrame2d.tcl
source EigenFrame.tcl
source EigenFrame.Extra.tcl
source AISC25.tcl

# Shells
source Shell/PinchedCylinder.tcl
# source Shell/PlanarShearWall.tcl

