
foreach model {
    Ex9a.build.UniaxialSection2D.tcl
    Ex9b.build.WSection2D.tcl
    Ex9c.build.RCSection.RectUnconfinedSymm2D.tcl
    Ex9d.build.RCSection.RectConfinedSymm2D.tcl
    Ex9e.build.RCSection.Rect2D.tcl
    Ex9f.build.RCSection.Circ2D.tcl
  } {
    wipe
    puts "Model: $model"
    source $model
    source Ex9.analyze.MomentCurvature2D.tcl
}

foreach model {
    Ex9a.build.UniaxialSection3D.tcl
    Ex9b.build.WSection3D.tcl
    Ex9c.build.RCSection.RectUnconfinedSymm3D.tcl
    Ex9d.build.RCSection.RectConfinedSymm3D.tcl
    Ex9e.build.RCSection.Rect3D.tcl
    Ex9f.build.RCSection.Circ3D.tcl
    Ex9g.build.HollowSection3D.tcl
  } {
    wipe
    puts "Model: $model"
    source $model
    source Ex9.analyze.MomentCurvature3D.tcl
}
