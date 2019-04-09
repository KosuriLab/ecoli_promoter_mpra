#!/bin/bash

# This script concatenates all qseq files with the same prefix into one large file and converts it to
# FASTQ format, only keeping those reads that pass filter. This script takes an optional third paramter, 
# '-c' to indicate if the qseq files are compressed
data=$1
prefix=$2

if [[ $3 = "-c" ]]; then
	echo "Processing $prefix ..."
	zcat $data/qseqs/$prefix*qseq.txt.gz | python qseq2fastq_stdin.py > $data/$prefix.fastq
else
	cat $data/qseqs/$prefix*qseq.txt | python qseq2fastq_stdin.py > $data/$prefix.fastq
fi


