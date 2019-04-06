prefix=$1
output=$2

echo "Aligning sample $prefix"
pear -j 30 --max-uncalled-base 0 -f $prefix*R2* -r $prefix*R1* -o ${prefix}.paired

mv ${prefix}.paired* $output
