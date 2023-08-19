#!/usr/bin/bash
#
build="/home/claudio/opensees/OpenSeesRT/build"
lib="/src/libg3/SRC/runtime/libOpenSeesRT.so"

opensees="python -m opensees"
#opensees="$build/temp.linux-x86_64-cpython-39_stack/src/libg3/SRC/executable/OpenSees"
perf="perf stat -B -r 5 -e cache-references,cache-misses,cycles,instructions,faults,page-faults,branches,branch-misses"

script="Ex8.genericFrame3D.tcl"

export OPENSEESRT_LIB="$build/temp.linux-x86_64-cpython-39_stack/$lib"
$perf $opensees $script

export OPENSEESRT_LIB="$build/temp.linux-x86_64-cpython-39_local/$lib"
$perf $opensees  $script

