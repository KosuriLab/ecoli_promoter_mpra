library(dplyr)
library(tidyr)
library(purrr)
library(Biostrings)

options(stringsAsFactors = F)
options(scipen = 10000)

args = commandArgs(trailingOnly=TRUE)

# args <- c('../../processed_data/frag/lb',
#           '../../processed_data/frag/barcode_maps.txt',
#           '../../processed_data/frag/frag_stats.txt',
#           '../../processed_data/frag/U00096.2_frag-rLP5_LB_expression.txt')


count_folder <- args[1]
bc_map_file <- args[2]
frag_ref_file <- args[3]
output_name <- args[4]

print(paste("Count folder:", count_folder))
print(paste("Barcode map file:", bc_map_file))
print(paste("Fragment reference file:", frag_ref_file))
print(paste("Output name:", output_name))

# Read in all Sequencing data for rLP5-Frag

# read in barcode count files
filelist = list.files(path = count_folder,
                      pattern = '^counts_*',
                      full.names = T)

sample_names <- c('DNA1_1', 'DNA1_2', 'DNA2_1', 'DNA2_2', 'RNA1_1', 'RNA1_2', 'RNA2_1', 'RNA2_2')

for(i in seq(1:length(filelist))) {
    sample_name <- sample_names[i]
    x <- read.table(filelist[i], col.names=c(sample_name, 'barcode'), header = F)
    x[[sample_name]] <- 1000000*x[[sample_name]]/sum(x[[sample_name]])  #Normalizes by RPM
    assign(sample_name, x)  
}

count_list <- list(DNA1_1, DNA1_2, DNA2_1, DNA2_2, RNA1_1, RNA1_2, RNA2_1, RNA2_2)

# combine reads for all barcodes 
# feed list into reduce to full join entire list of data frames
all_counts <- count_list %>% purrr::reduce(full_join, by = "barcode")
print(paste("Number of unique sequenced barcodes, all replicates in DNA and RNA:", nrow(all_counts)))

# remove individual sample files
rm(count_list, x)
rm(list = sample_names)

all_counts[is.na(all_counts)] <- 0

mapped_fragstats <- read.table(frag_ref_file, header = F, fill = T, 
                               col.names = c('barcode', 'start', 'end', 'strand', 'fragment'))
mapped_fragstats <- mapped_fragstats[!duplicated(mapped_fragstats$barcode),]

fragstats <- left_join(mapped_fragstats, all_counts, by = 'barcode')

frag <- fragstats %>%
    group_by(fragment) %>% 
    mutate(num_mapped_barcodes = n()) %>%
    filter(DNA1_1 > 0 | DNA1_2 > 0 | DNA2_1 > 0 | DNA2_2 > 0) %>%
    mutate(num_integrated_barcodes = n()) %>%
    filter(num_integrated_barcodes >= 1) %>%
    mutate(RNA_exp_1 = (sum(RNA1_1)+sum(RNA1_2))/(sum(DNA1_1)+sum(DNA1_2)),
           RNA_exp_2 = (sum(RNA2_1)+sum(RNA2_2))/(sum(DNA2_1)+sum(DNA2_2)),
           RNA_exp_ave = ((RNA_exp_1 + RNA_exp_2)/2),
           DNA_sum_1 = sum(DNA1_1)+sum(DNA1_2),
           DNA_sum_2 = sum(DNA1_1)+sum(DNA1_2),
           DNA_ave = ((DNA_sum_1 + DNA_sum_2)/2)) %>% 
    ungroup() %>% 
    # Filter out fragments with low expression
    filter(RNA_exp_1 >= .1 & RNA_exp_2 >= .1, RNA_exp_1 < 10000 & RNA_exp_2 < 10000, DNA_ave>1) %>% 
    mutate(variation = abs(log2(RNA_exp_2/RNA_exp_1))) %>% 
    # Remove fragments that have a 5x difference between biological replicates
    filter(variation < 2.32) %>% 
    select(fragment, RNA_exp_1, RNA_exp_2, RNA_exp_ave, DNA_sum_1, DNA_sum_2, DNA_ave, num_mapped_barcodes, num_integrated_barcodes, start, end, strand, variation) %>% 
    distinct() 

write.table(frag, file = output_name, quote = F, row.names = F)

