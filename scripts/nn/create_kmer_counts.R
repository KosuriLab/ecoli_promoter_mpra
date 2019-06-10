load_pkgs <- function(pkgs){
    new_pkgs <- pkgs[!(pkgs %in% installed.packages()[, 'Package'])]
    if(length(new_pkgs)) install.packages(new_pkgs)
    for(pkg in pkgs){
        suppressWarnings(suppressMessages(library(pkg, character.only = T)))
    }
}

pkgs <- c('kmer', 'dplyr', 'seqinr')
load_pkgs(pkgs)


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

