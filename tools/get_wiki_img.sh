

while read line; do
 curl -L "https://opensees.berkeley.edu/wiki/index.php/Special:FilePath/$line" > "$1/$line"
done


