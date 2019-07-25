library(dplyr)

Endo2 <- read.table("../../processed_data/endo_tss/lb/rLP5_Endo2_lb_expression.txt", 
                    header = T)

# identify two standard deviations greater than median of negatives
neg <- subset(Endo2, grepl("neg_control", Endo2$name))
neg_median <- median(neg$RNA_exp_ave)
neg_sd <- sd(neg$RNA_exp_ave)

# Subset all promoters that are greater than 3sd from the mean
positive_Endo2 <- filter(Endo2, RNA_exp_ave > (neg_median+3*(neg_sd)))
negative_Endo2 <- filter(Endo2, RNA_exp_ave < (neg_median+3*(neg_sd)))

# Write out active and inactive promoters
positive_Endo2 %>%
    mutate(chrom = 'U00096.2') %>%
    select(chrom, start, end, name, RNA_exp_ave, strand) %>%
    filter(start > 0) %>%
    write.table('../../processed_data/endo_tss/lb/tss_positives.bed', 
                sep = '\t', col.names = F, row.names = F, quote = F)

negative_Endo2 %>%
    mutate(chrom = 'U00096.2') %>%
    select(chrom, start, end, name, RNA_exp_ave, strand) %>%
    filter(start > 0) %>%
    write.table('../../processed_data/endo_tss/lb/tss_negatives.bed', 
                sep = '\t', col.names = F, row.names = F, quote = F)

# convert to FASTA
system(paste('bedtools getfasta',
             '-fi ../../ref/U00096.2.fasta',
             '-bed ../../processed_data/endo_tss/lb/tss_positives.bed',
             '-fo ../../processed_data/endo_tss/lb/tss_positives.fasta',
             '-name -s'))

system(paste('bedtools getfasta',
             '-fi ../../ref/U00096.2.fasta',
             '-bed ../../processed_data/endo_tss/lb/tss_negatives.bed',
             '-fo ../../processed_data/endo_tss/lb/tss_negatives.fasta',
             '-name -s'))
