DIR=$1
echo "Extracting & counting barcodes..."

for i in ${DIR}/*.fastq.gz; 
do 
	zcat $i |
	awk 'NR%4==2' |
	cut -c 1-20 | 
	rev | 
	tr ACGT TGCA |
	sort -T ./ --parallel=20 | 
	uniq -c > counts_$(basename $i .fastq.gz).txt; 
done
