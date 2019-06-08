library(kmer)
library(dplyr)
library(seqinr)
options(stringsAsFactors = F)

args = commandArgs(trailingOnly=TRUE)

filename <- args[1]
outfile <- args[2]

print(filename)

# read in data, assume first column is sequence

data <- read.table(filename)
sequence <- data[,1]

seq_list <- strsplit(sequence, '')
counts <- as.data.frame(kcount(seq_list, k = 6))

output <- bind_cols(data, counts)

write.table(output, outfile, row.names = F, quote = F, sep = '\t', col.names = F)

