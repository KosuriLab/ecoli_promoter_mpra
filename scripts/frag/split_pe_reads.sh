prefix=$1

echo "Combining reads from sample ${prefix}..."

awk  '{gsub("GCATGTGAGACCGGATGCTAACTAAACACCGCTAGC","\t",$0); print;}' $prefix*R1* | \
paste $prefix*R2* - | \
awk 'NR%4==2' | \
awk '{print $1, $3, $2}' | \
awk 'NF==3{print}{}' > $(basename $prefix).txt
