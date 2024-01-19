
foreach model {
    ElasticElement
    InelasticFiberSection
    InelasticSection
  } {
  puts "Model: $model"
  foreach analysis {
    Dynamic.EQ.bidirect
    Dynamic.EQ.multipleSupport
    Dynamic.EQ.Uniform
    Dynamic.sine.multipleSupport
    Dynamic.sine.Uniform
    Static.Cycle
    Static.Push
  } {
    puts "       $analysis"
    source Ex4.Portal2D.build.$model.tcl
    source Ex4.Portal2D.analyze.$analysis.tcl
  }
}
