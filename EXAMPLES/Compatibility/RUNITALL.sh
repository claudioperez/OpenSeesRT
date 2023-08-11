#!/bin/bash
#
Results="README.md"
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

