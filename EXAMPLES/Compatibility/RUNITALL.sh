#!/bin/bash
#
Results="README.md"

export OPENSEESRT_LIB="/home/claudio/opensees/OpenSeesRT/build/temp.linux-x86_64-cpython-39_stack/src/libg3/SRC/runtime/libOpenSeesRT.so"

OPENSEES="python -m opensees"
#OPENSEES=OpenSees

cat - > $Results <<EOF
# Examples

| Status  |     File     |
|---------|--------------|
EOF

# -e Shell_Dynamic
#
ls $@ | grep -v -e "Rocking" -e "Brick27"  | while read i; do
  echo "FILE $i"; 
  $OPENSEES $i; 
  printf "| $?\t|"' `'"$i"'` |\n' | tee -a $Results;
done

