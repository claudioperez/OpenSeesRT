set files {
      TestFPS2d_0.tcl
      TestFPS2d_1.tcl
      TestFPS2d_2.tcl
      TestFPS2d_3.tcl
      TestFPS2d_4.tcl
      TestFPS3d_0.tcl
      TestFPS3d_1.tcl
      TestFPS3d_2.tcl
      TestFPS3d_3.tcl
      TestFPS3d_4.tcl
      TestSlider2d_0.tcl
      TestSlider2d_1.tcl
      TestSlider2d_2.tcl
      TestSlider2d_3.tcl
      TestSlider2d_4.tcl
      TestSlider3d_0.tcl
      TestSlider3d_1.tcl
      TestSlider3d_2.tcl
      TestSlider3d_3.tcl
      TestSlider3d_4.tcl
}

foreach file $files {
  wipe
  puts "$file\t[source $file]"
}
