
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
Ex5.Frame2D.analyze.Dynamic.EQ.bidirect.tcl
Ex5.Frame2D.analyze.Dynamic.EQ.multipleSupport.tcl
Ex5.Frame2D.analyze.Dynamic.EQ.Uniform.tcl
Ex5.Frame2D.analyze.Dynamic.sine.multipleSupport.tcl
Ex5.Frame2D.analyze.Dynamic.sine.Uniform.tcl
Ex5.Frame2D.analyze.Static.Cycle.tcl
Ex5.Frame2D.analyze.Static.Push.tcl
Ex5.Frame2D.build.ElasticSection.tcl
Ex5.Frame2D.build.InelasticFiberRCSection.tcl
Ex5.Frame2D.build.InelasticFiberWSection.tcl
Ex5.Frame2D.build.InelasticSection.tcl
test.tcl
