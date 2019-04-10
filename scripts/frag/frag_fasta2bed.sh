DIR=$1
# Get fasta sequence from bed files, need to do + strand maps and then minus strand
echo "Extracting full sequences and barcode..."

fasta=../../ref/Escherichia_coli_K-12_MG1655.fasta

awk '$3==99' $DIR/sam_stats.txt | 
awk '{print "U00096.2", $4-1, $4+$6-1, $8, "1", "+"}' | 
tail -n +2 | 
sed 's/  */\t/g' > Gfrag_plus.bed

bedtools getfasta -s -fi $fasta -bed Gfrag_plus.bed -fo FragmentBarcode_plus.fa.out -name -tab

## Minus Strand ##
awk '$3==83' $DIR/sam_stats.txt | 
awk '{print "U00096.2", $5-1, $5-1-$6, $8, "1", "-"}' | 
tail -n +2 | 
sed 's/  */\t/g' > Gfrag_minus.bed

bedtools getfasta -s -fi $fasta -bed Gfrag_minus.bed -fo FragmentBarcode_minus.fa.out -name -tab

cat FragmentBarcode_plus.fa.out FragmentBarcode_minus.fa.out | 
sort --parallel=30 -k 1 | 
uniq > FragmentBarcodes.txt

cat Gfrag_plus.bed Gfrag_minus.bed | 
sort --parallel=30 -k 4 | 
uniq | 
join -1 4 - FragmentBarcodes.txt | 
awk '{print $1, $3, $4, $6, $7}' > $DIR/frag_stats.txt

rm -f FragmentBarcode_*.fa.out
rm -f Gfrag*
echo "Done!"
