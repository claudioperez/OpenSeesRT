
foreach model {
  ElasticElement 
  InelasticFiberSection
  InelasticSection
  } {

  puts "Model: $model"
  foreach analysis {
    Static.Push
    Dynamic.EQ.Uniform
    } {
    puts "       $analysis"
    wipe
    source Ex3.Canti2D.build.$model.tcl
    source Ex3.Canti2D.analyze.$analysis.tcl
  }
}
