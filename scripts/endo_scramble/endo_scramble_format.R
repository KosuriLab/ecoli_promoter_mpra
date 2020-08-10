library(tidyverse)
#library(tidyr)

options(stringsAsFactors = F)

args = commandArgs(trailingOnly=TRUE)

# infile <- args[1]
# outfile <- args[2]
infile <- '../../processed_data/endo_scramble/endo_scramble_expression.txt'
outfile <- '../../processed_data/endo_scramble/endo_scramble_expression_formatted.txt'

data <- read.table(file = infile, header = T)


data <- data %>% 
    mutate(name = gsub('>', '', orig_name),
           name = gsub('_rc', '', name)) %>% 
    # parse name
    separate(name, into = c('tss_name', 'tss_position', 'strand_scramble_loc'), sep = ',',
             convert = T, remove = F) %>% 
    separate(strand_scramble_loc, into = c('strand', 'scramble_loc'), sep = '_') %>% 
    mutate(scramble_loc = gsub('unscrambled', NA, scramble_loc),
           scramble_loc = gsub('scrambled', '', scramble_loc)) %>% 
    separate(scramble_loc, into = c('scramble_start', 'scramble_end'),
             sep = '-', convert = T) %>% 
    # calculate genomic position of scramble start and end
    mutate(var_left = ifelse(strand == '+', tss_position - 120, tss_position - 30),
           var_right = ifelse(strand == '+', tss_position + 30, tss_position + 120),
           scramble_start_pos = ifelse(strand == '+',
                                       var_left + scramble_start, 
                                       var_right - scramble_start),
           scramble_end_pos = ifelse(strand == '+', 
                                     var_left + scramble_end, 
                                     var_right - scramble_end)) %>%
           #scramble_start_pos = var_left + scramble_start,
           #scramble_end_pos = var_left + scramble_end) %>%
    # calculate scramble position relative to TSS
    mutate(scramble_pos_rel_tss = ifelse(is.na(scramble_start_pos),
                                         NA,
                                         ifelse(strand == '+', 
                                                scramble_start_pos - tss_position,
                                                tss_position - scramble_start_pos))) %>% 
    select(name:scramble_pos_rel_tss, variant, RNA_exp_sum_1_1:num_barcodes_integrated)


data$category <- "scramble"
data$category[grep("pos_control", data$name)] <- "pos_control"
data$category[grep("neg_control", data$name)] <- "neg_control"
data$category[grep("unscrambled", data$name)] <- "unscrambled"

write.table(data, file = outfile, quote = F, row.names = F)
