library(kmer)
library(dplyr)
library(seqinr)
options(stringsAsFactors = F)

args = commandArgs(trailingOnly=TRUE)

filename <- args[1]
outfile <- args[2]

print(filename)

# filename <- '../../processed_data/combined/tss_scramble_peak_expression_model_format_test_genome_split_negatives.fasta'
# filename <- '../../processed_data/combined/tss_scramble_peak_expression_model_format_floored_train_genome_split.txt'

# read in data, assume first column is sequence
if(endsWith(filename, 'txt')){
    data <- read.table(filename)
    sequence <- data[,1]
} else if(endsWith(filename, 'fasta')) {
    sequence <- toupper(unlist(read.fasta(filename, as.string = T)))
} else {
    stop("Please specify .txt or .fasta")
}

seq_list <- strsplit(sequence, '')
counts <- as.data.frame(kcount(seq_list, k = 6))

if(endsWith(filename, 'txt')) {
    output <- bind_cols(data, counts)
}
if (endsWith(filename, 'fasta')) {
    output <- data.frame(sequence, counts)
}

write.table(output, outfile, row.names = F, quote = F, sep = '\t', col.names = F)

