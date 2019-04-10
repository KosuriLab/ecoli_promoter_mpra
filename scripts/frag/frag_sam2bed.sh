DIR=$1

# Turn mapping into bed files with barcodes as names
echo "Converting mapping into bed files..."
sam=$DIR/U00096.2_Align.sam
awk '{print $1,$10,$2,$4,$8,$9}' $sam | awk 'NR%2==0' | tail -n +2 > Map5p.txt
awk '{print $1,$10}' $sam | awk 'NR%2==1' | tail -n +3 > Map3p.txt
join -a1 -e 0 Map5p.txt Map3p.txt | sort -n -k1 | paste - $DIR/barcode_maps.txt > $DIR/sam_stats.txt
rm -f Map5p.txt Map3p.txt