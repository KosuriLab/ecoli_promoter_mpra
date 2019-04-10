#!/bin/bash

DIR=$1
OUTPUT=$2

echo "Extracting & counting barcodes..."

for i in ${DIR}/rLP5*.fastq.gz; 
do 
	awk 'NR%4==2' $i | 
	cut -c 1-20 |
	sort --parallel=30 |
	uniq -c > ${OUTPUT}/$(basename $i .fastq.gz).txt; 
done

# pretty names
for file in ${OUTPUT}/rLP5*.txt;
do
	mv $file ${file/_S*/}.txt;
done
