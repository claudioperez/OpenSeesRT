#!/bin/bash
#
Results="README.md"

cat - > $Results <<EOF
# Examples

| Status  |     File     |  Message              |
|---------|--------------|-----------------------|
EOF

ls $@ | grep -v -e "Proc.*" -e '.*FrictionPendulum.*' -e lasticPileSection.tcl -e Dispwall1-cg.tcl -e Dynamic.EQ.Uniform_LimitState.tcl -e Ex6.genericFrame[23]D.build.InelasticFiber -e "Lib.*" -e TestSlider -e "Read.*" -e "Ex.*\.analyze\..*" -e "^PUL.*" -e "CenterCol[A-Z]" -e '.*_input_.*'  -e Tags.tcl -e 'NRHA_IR.tcl' -e 'NR94cnp.tcl' | while read i; do
  echo "FILE $i"; 
  msg="$(python -m opensees $i 2>&1 >/dev/null)";
  code="$?"
  echo "$msg"
  msg=$(tail -1 <<< "$msg")
  printf "|  $code  |"' `'"$i"'`'" | $msg |\n"  >> $Results; #| tee -a $Results;
done

