if 1 {
  puts " -------------Static Pushover Analysis -------------"
  foreach sec {RCsec Wsec} {
    wipe
    source Ex8.genericFrame3D.build.$sec.tcl
    source Ex8.genericFrame3D.analyze.Static.Push.tcl
  }
}
# To run RC Model, Uniform Earthquake Excitation
if 0 {
  wipe
  puts " -------------Uniform Earthquake Excitation -------------"
  source Ex8.genericFrame3D.build.RCsec.tcl
  source Ex8.genericFrame3D.analyze.Dynamic.EQ.Uniform.tcl
}
