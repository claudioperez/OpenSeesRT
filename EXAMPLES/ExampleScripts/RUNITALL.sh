
Results="README.md"

cat - > $Results <<EOF
# Examples

| Status  |     File     |
|---------|--------------|
EOF

ls $@ | grep -v -e "Rocking" -e Brick27 | while read i; do
  echo "FILE $i"; 
  python -m opensees $i; 
  printf "| $?\t|"' `'"$i"'` |\n' | tee -a $Results;
done

