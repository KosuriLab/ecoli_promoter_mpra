#get only pairs of reads over 1bp and trim A-tails, change print and name for 3prime end
echo "Extracting reads and barcode..."

DATA=$1

cat $DATA/pFrag_80.txt $DATA/pFrag_83.txt | \
awk 'NF==3' | \
awk 'length($1)>1' | \
awk 'length($2)>1' | \
awk '{print $1}' | \
cut -c 2- | \
awk 'NF==1' > $DATA/5prime.txt 


cat $DATA/pFrag_80.txt $DATA/pFrag_83.txt | \
awk 'NF==3' | \
awk 'length($1)>1' | \
awk 'length($2)>1' | \
awk '{print $2}' | \
cut -c 2- | \
awk 'NF==1' > $DATA/3prime.txt

cat $DATA/pFrag_80.txt $DATA/pFrag_83.txt | \
awk 'NF==3' | \
awk 'length($1)>1' | \
awk 'length($2)>1' | \
awk '{print $1,$3}' | \
cut -c 2- | \
awk '{print $2, $1}' > $DATA/barcode_maps.txt
