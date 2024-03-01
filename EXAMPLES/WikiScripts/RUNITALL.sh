#!/bin/bash
#
Results="Status.md"
# Results="Memory.md"
# export LD_PRELOAD=/home/claudio/mambaforge/envs/py39/x86_64-conda-linux-gnu/lib/libasan.so
export OPENSEESRT_LIB=/home/claudio/packages/opensees-pypi/build/temp.linux-x86_64-cpython-39_debug/src/libg3/SRC/runtime/libOpenSeesRT.so
# export OPENSEESRT_LIB="/home/claudio/opensees/OpenSeesRT/build/temp.linux-x86_64-cpython-39_debug/src/libg3/SRC/runtime/libOpenSeesRT.so"

cat - > $Results <<EOF
# Examples

| Status  |     File     |  Message              |
|---------|--------------|-----------------------|
EOF

for i in */test.tcl; do
  dir="$(dirname $i)"
  cd $dir
  echo "FILE $i"; 
  msg="$(python -m opensees test.tcl 2>&1 >/dev/null)";
  code="$?"
  echo "$msg"
  msg=$(tail -1 <<< "$msg")
  printf "|  $code  | "'[`'"$dir"'`]('./$dir") | $msg |\n"  >> ../$Results; #| tee -a $Results;
  #
  cd ../
done
#exit
ls *.tcl | grep -v -e Ex8 -e Dynamic.EQ.Uniform_LimitState -e Ex[68].genericFrame[23]D.build'\.' -e "Lib.*" -e "Read.*" -e "Ex.*\.analyze\..*" -e "CenterCol[A-Z]" -e '.*_input_.*'  -e Tags.tcl  | while read i; do
  echo "FILE $i"; 
# msg="$(OpenSees $i 2>&1 >/dev/null)";
  msg="$(python -m opensees $i 2>&1 >/dev/null)";
# msg="$(python -m opensees $i)";
  code="$?"
  echo "$msg"
  msg=$(tail -1 <<< "$msg")
# printf "|  $code  |"' `'"$i"'`'" | $msg |\n"  >> $Results; #| tee -a $Results;
  printf "|  $code  | "'[`'"$i"'`]('./$i") | $msg |\n"  >> $Results;
done

