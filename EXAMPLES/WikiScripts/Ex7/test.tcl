
foreach section {RCsec Wsec} {
  foreach analysis {Dynamic.EQ.bidirect
                    Dynamic.EQ.Uniform
                    Dynamic.EQ.multipleSupport
                    Static.Cycle
                    Static.Push
                  } {
#                   Dynamic.sine.multipleSupport
    puts "\t$section \t $analysis"
    source Ex7.Frame3D.build.$section.tcl
    source Ex7.Frame3D.analyze.$analysis.tcl

  }
}
