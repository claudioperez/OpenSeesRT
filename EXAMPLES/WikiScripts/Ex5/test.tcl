
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
                    Dynamic.sine.Uniform
                    Static.Cycle
                    Static.Push
    } {
    puts "       $analysis"
    wipe
    source Ex5.Frame2D.build.$model.tcl
    source Ex5.Frame2D.analyze.$analysis.tcl
  }
}
