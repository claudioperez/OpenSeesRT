
foreach model {
    ElasticSection
    InelasticFiberRCSection
    InelasticFiberWSection
    InelasticSection
  } {

  puts "Model: $model"
  foreach analysis {
    Dynamic.EQ.bidirect
    Dynamic.EQ.multipleSupport
    Dynamic.EQ.Uniform
    Dynamic.sine.multipleSupport
    Static.Cycle
    Static.Push
    } {
    puts "       $analysis"
    wipe
    source Ex6.genericFrame2D.build.$model.tcl
    source Ex6.genericFrame2D.analyze.$analysis.tcl
  }
}
