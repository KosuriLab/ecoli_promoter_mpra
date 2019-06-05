# source('https://bioconductor.org/biocLite.R')
# biocLite('GenomicRanges')
# biocLite('rtracklayer')
# biocLite('BSgenome')
# install.packages('ROCR')
# install.packages('kernlab')
# install.packages('seqinr')
# install.packages('gkmSVM')

library('gkmSVM')
library('seqinr')
library('ggplot2')
library('dplyr')
require(cowplot)

options(stringsAsFactors = F)
set.seed(123)


# run gkmSVM on TSS library sequences
train <- read.table('../../processed_data/endo_tss/lb/model_files/tss_expression_model_format_train_genome_split.txt',
                    header = F, sep = '\t', col.names = c('sequence', 'expression'))
train$name <- paste0('seq', seq(1:nrow(train)))

# split train into positive and negative
train_positive <- filter(train, expression >= 1)
train_negative <- filter(train, expression < 1)

write.fasta(as.list(train_positive$sequence), names = train_positive$name, nbchar = 300,
            file.out = '../../processed_data/endo_tss/lb/model_files/tss_train_genome_split_positives.fasta')
write.fasta(as.list(train_negative$sequence), names = train_negative$name, nbchar = 300,
            file.out = '../../processed_data/endo_tss/lb/model_files/tss_train_genome_split_negatives.fasta')

kmer_length = 10
ungapped <- 8
# print("Creating kernel...")
# gkmsvm_kernel(posfile = "../../processed_data/endo_tss/lb/model_files/tss_train_genome_split_positives.fasta", 
#               negfile = "../../processed_data/endo_tss/lb/model_files/tss_train_genome_split_negatives.fasta", 
#               outfile = '../../processed_data/endo_tss/lb/model_files/gkmsvm_10mer_8ungapped_kernel', 
#               L = kmer_length, K = ungapped)
# 
# print("Training...")
# # perform SVM training with cross-validation
# result <- gkmsvm_trainCV(kernelfn = '../../processed_data/endo_tss/lb/model_files/gkmsvm_10mer_8ungapped_kernel',
#                          posfn= "../../processed_data/endo_tss/lb/model_files/tss_train_genome_split_positives.fasta", 
#                          negfn = "../../processed_data/endo_tss/lb/model_files/tss_train_genome_split_negatives.fasta",
#                          svmfnprfx='../../processed_data/endo_tss/lb/model_files/gkmsvm_10mer_8ungapped_kernel', 
#                          outputCVpredfn='../../processed_data/endo_tss/lb/model_files/gkmsvm_10mer_8ungapped_cvpred',
#                          outputROCfn='../../processed_data/endo_tss/lb/model_files/gkmsvm_10mer_8ungapped_roc', 
#                          L = kmer_length, K = ungapped, showPlots = T, 
#                          outputPDFfn = '../../processed_data/endo_tss/lb/model_files/gkmsvm_10mer_8ungapped_curves')
# 
# ggsave('../../processed_data/endo_tss/lb/model_files/gkmsvm_s10mer_8ungapped_ROC_PR_curves.png')

print("Classifying test set...")
# read in test set and convert to FASTA
test <- read.table('../../processed_data/endo_tss/lb/model_files/tss_expression_model_format_test_genome_split.txt',
                   header = F, sep = '\t', col.names = c('sequence', 'expression'))
test$name <- paste0('seq', seq(1:nrow(test)))
write.fasta(as.list(test$sequence), names = test$name, nbchar=300,
            file.out = '../../processed_data/endo_tss/lb/model_files/tss_expression_model_format_test_genome_split.fasta')

gkmsvm_classify(seqfile = '../../processed_data/endo_tss/lb/model_files/tss_expression_model_format_test_genome_split.fasta',
                svmfnprfx = '../../processed_data/endo_tss/lb/model_files/gkmsvm_10mer_8ungapped_kernel',
                outfile = '../../processed_data/endo_tss/lb/model_files/gkmsvm_10mer_8ungapped_kernel_test_results',
                L = kmer_length, K = ungapped)


