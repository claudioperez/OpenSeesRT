    puts " -------------Elastic Model -------------"
    puts " -------------Static Pushover Analysis -------------"
    source Ex7.Frame3D.build.Wsec.tcl
    source Ex7.Frame3D.analyze.Static.Push.tcl

#   To run RC Model, Uniform Earthquake Excitation

    puts " -------------Uniaxial Inelastic Section, Nonlinear Model -------------"
    puts " -------------Uniform Earthquake Excitation -------------"
    source Ex7.Frame3D.build.RCsec.tcl
    source Ex7.Frame3D.analyze.Dynamic.EQ.Uniform.tcl

