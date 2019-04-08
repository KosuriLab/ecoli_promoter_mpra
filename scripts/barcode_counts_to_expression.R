library(dplyr)
library(tidyr)
library(purrr)
library(Biostrings)

options(stringsAsFactors = F)
options(scipen = 10000)

args = commandArgs(trailingOnly=TRUE)

# args <- c('../../processed_data/endo_tss/lb',
#           '../../processed_data/endo_tss/lb/endo_lb_bc_map.txt',
#           '../../processed_data/endo_tss/lb/endo_lb_controls_bc_map.txt',
#           '../../ref/endo_lib_2016_controls_clean.txt',
#           '../../processed_data/endo_tss/lb/rLP5_Endo2_lb_expression.txt')

# args <- c('../../processed_data/endo_scramble',
#           '../../processed_data/endo_scramble/endo_scramble_combined_bc_map.txt',
#           '../../processed_data/endo_scramble/endo_scramble_combined_controls_bc_map.txt',
#           '../../ref/20180507_active_tss_scrambled10_stride5.txt',
#           '../../processed_data/endo_scramble/endo_scramble_expression.txt')
count_folder <- args[1]
bc_map_file <- args[2]
bc_map_file_controls <- args[3]
ref_file <- args[4]
output_name <- args[5]

print(paste("Count folder:", count_folder))
print(paste("Barcode map file:", bc_map_file))
print(paste("Barcode map control file", bc_map_file_controls))
print(paste("Reference file:", ref_file))
print(paste("Output name:", output_name))

# read in barcode count files
filelist = list.files(path = count_folder,
                      pattern = '^counts_*',
                      full.names = T)
# should be named so files are in following order, 1st number is biological replicate
# and second number is technical replicate
if(length(filelist) == 8) {
    sample_names <- c('DNA1_1', 'DNA1_2', 'DNA2_1', 'DNA2_2', 'RNA1_1', 'RNA1_2', 'RNA2_1', 'RNA2_2')
} else if(length(filelist) == 4) {
    sample_names = c('DNA1', 'DNA2', 'RNA1', 'RNA2')
} else {
    stop("Unrecognized experiment design, do not know how to process")
}

for(i in seq(1:length(filelist))) {
    sample_name <- sample_names[i]
    x <- read.table(filelist[i], col.names=c(sample_name, 'barcode'), header = F)
    x[[sample_name]] <- 1000000*x[[sample_name]]/sum(x[[sample_name]])  #Normalizes by RPM
    assign(sample_name, x)  
}

if(length(filelist) == 8) {
    count_list <- list(DNA1_1, DNA1_2, DNA2_1, DNA2_2, RNA1_1, RNA1_2, RNA2_1, RNA2_2)
} else if(length(filelist) == 4) {
    count_list <- list(DNA1, DNA2, RNA1, RNA2)
} 

# combine reads for all barcodes 
# feed list into reduce to full join entire list of data frames
all_counts <- count_list %>% purrr::reduce(full_join, by = "barcode")
print(paste("Number of unique sequenced barcodes, all replicates in DNA and RNA:", nrow(all_counts)))

# # rearrange columns
# all_counts <- select(all_counts, barcode, DNA1_1, DNA1_2:RNA2_2)

# remove individual sample files
rm(count_list, x)
rm(list = sample_names)

# read in mapping files
bc_map <- read.table(bc_map_file, header = F, sep = ',', 
                       col.names = c('barcode', 'variant', 'count'))

# read in reference so we can label mapping file
ref <- read.table(ref_file, col.names = c('name', 'sequence'))

# add reverse complement for each sequence so we can properly label everything
ref_with_rc <- ref %>% 
    mutate(name = paste0(name, '_rc'),
           sequence = as.character(reverseComplement(DNAStringSet(sequence)))) %>% 
    select(name, sequence) %>% 
    bind_rows(select(ref, name, sequence)) %>% 
    # trim primers so sequence is 150 to match sequence in barcode map
    mutate(variant = toupper(substr(sequence, 25, 174)))

# add variant name from reference
bc_map <- bc_map %>% 
    left_join(select(ref_with_rc, -sequence), by = 'variant')

# add in positive controls
bc_map_controls <- read.table(bc_map_file_controls, header = T, 
                           col.names = c('barcode', 'num_unique_variant', 'count',
                                         'count_most_common', 'variant', 'name'))
bc_map <- bind_rows(bc_map, select(bc_map_controls, barcode, variant, count, name))
print(paste("Number of barcodes in mapping run:", nrow(bc_map)))

# only keep those that match reference
bc_map <- filter(bc_map, !is.na(name))
print(paste("Number of barcodes in mapping run and mapped to reference:", nrow(bc_map)))

# only keep barcode counts for mapped barcodes
mapped_bc_counts <- inner_join(bc_map, all_counts, by = 'barcode')
# replace NA with 0
mapped_bc_counts[is.na(mapped_bc_counts)] <- 0
print(paste("Number of mapped and sequenced barcodes, all replicates in DNA and RNA:", nrow(mapped_bc_counts)))

# compute summary statistics
if(length(filelist) == 8) {
    expression <- mapped_bc_counts %>% 
        group_by(variant) %>% 
        mutate(num_barcodes_mapped = n()) %>%
        filter(DNA1_1 > 0, DNA1_2 > 0, DNA2_1 > 0, DNA2_2 > 0) %>%
        mutate(num_barcodes_integrated = n()) %>%
        # filter out promoters with fewer than 3 barcodes integrated
        filter(num_barcodes_integrated >= 3) %>% 
        mutate(RNA_exp_sum_1_1 = sum(RNA1_1)/(sum(DNA1_1)),
               RNA_exp_sum_1_2 = sum(RNA1_2)/(sum(DNA1_2)),
               RNA_exp_sum_2_1 = sum(RNA2_1)/(sum(DNA2_1)),
               RNA_exp_sum_2_2 = sum(RNA2_2)/(sum(DNA2_2)),
               RNA_exp_sum_1 = ((RNA_exp_sum_1_1+RNA_exp_sum_1_2)/2),
               RNA_exp_sum_2 = ((RNA_exp_sum_2_1+RNA_exp_sum_2_2)/2),
               RNA_exp_sum_ave = ((RNA_exp_sum_1+RNA_exp_sum_2)/2),
               DNA_1 = (sum(DNA1_1)+sum(DNA1_2)),
               DNA_2 = (sum(DNA2_1)+sum(DNA2_2)),
               DNA_ave = ((DNA_1 + DNA_2)/2),
               expn_1_1 = RNA1_1 / DNA1_1,
               expn_1_2 = RNA1_2 / DNA1_2,
               expn_2_1 = RNA2_1 / DNA2_1,
               expn_2_2 = RNA2_2 / DNA2_2,
               expn_med1 = median(c(expn_1_1, expn_1_2), na.rm = T),
               expn_med2 = median(c(expn_2_1, expn_2_2), na.rm = T),
               expn_med = (expn_med1 + expn_med2)/2) %>%
        ungroup() %>%
        select(variant, orig_name = name, RNA_exp_sum_1_1, RNA_exp_sum_1_2, RNA_exp_sum_2_1, RNA_exp_sum_2_2,
               RNA_exp_sum_1, RNA_exp_sum_2, RNA_exp_sum_ave, DNA_1, DNA_2, DNA_ave, expn_med,
               num_barcodes_mapped, num_barcodes_integrated) %>% 
        distinct() 
} else if(length(filelist) == 4) {
    expression <- mapped_bc_counts %>% 
        group_by(variant) %>% 
        mutate(num_barcodes_mapped = n()) %>%
        filter(DNA1 > 0, DNA2 > 0) %>%
        mutate(num_barcodes_integrated = n()) %>%
        # filter out promoters with fewer than 3 barcodes integrated
        filter(num_barcodes_integrated >= 3) %>% 
        mutate(DNA_sum = (sum(DNA2)+sum(DNA1)),
               RNA_exp_1 = sum(RNA1)/sum(DNA1),
               RNA_exp_2 = sum(RNA2)/sum(DNA2),
               RNA_exp_sum_ave = ((RNA_exp_1+RNA_exp_2)/2),
               expn_med1 = median(RNA1/DNA1, na.rm = T),
               expn_med2 = median(RNA2/DNA2, na.rm = T),
               expn_med = (expn_med1 + expn_med2) / 2) %>% 
        select(variant, orig_name = name, RNA_exp_1, RNA_exp_2, RNA_exp_sum_ave, DNA_sum, 
               expn_med1, expn_med2, expn_med, 
               num_barcodes_mapped, num_barcodes_integrated) %>% 
        distinct() 
}

percent_covered <- nrow(expression) / nrow(ref)
print(paste0("Mapped and integrated library: ", round(percent_covered,2) * 100, '%'))

# find threshold - 2 standard deviations greater than median negative control
negatives <- filter(expression, grepl('neg', orig_name))
neg_sd_sum <- sd(negatives$RNA_exp_sum_ave)
neg_median_sum <- median(negatives$RNA_exp_sum_ave, na.rm = T)
threshold_sum <- 2*neg_sd_sum + neg_median_sum
print(paste("Active threshold (sum):", threshold_sum))
print("Standardize so threshold = 1")

neg_sd_med <- sd(negatives$expn_med)
neg_median_med <- median(negatives$expn_med, na.rm = T)
threshold_med <- 2*neg_sd_med + neg_median_med
print(paste("Active threshold (median):", threshold_med))
print("Standardize so threshold = 1")

expression <- expression %>% 
    mutate(expn_sum_std = RNA_exp_sum_ave - threshold_sum,
           expn_med_std = expn_med - threshold_med)

write.table(expression, output_name, quote = F, row.names = F)


