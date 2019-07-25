args = commandArgs(trailingOnly = TRUE)
options(stringsAsFactors=FALSE)

filename = args[1]
output_filename = args[2]

oligos = read.table(file=filename, header=F, sep='\t')
colnames(oligos) = c('header', 'sequence')


# remove duplicate sequences
no_duplicates = oligos[!duplicated(oligos$sequence),]

write.table(no_duplicates, output_filename, quote=F, col.names=F, row.names=F, sep='\t')
