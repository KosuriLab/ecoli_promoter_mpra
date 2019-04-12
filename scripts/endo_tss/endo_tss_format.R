library(dplyr)
library(tidyr)

options(stringsAsFactors = F)

# args = commandArgs(trailingOnly=TRUE)
# 
# infile <- args[1]
# outfile <- args[2]
infile <- '../../processed_data/endo_tss/lb/rLP5_Endo2_lb_expression.txt'
outfile <- '../../processed_data/endo_tss/lb/rLP5_Endo2_lb_expression_formatted.txt'

# infile <- '../../processed_data/endo_tss/alt_landing_pads/fLP3/fLP3_Endo2_lb_expression.txt'
# outfile <- '../../processed_data/endo_tss/alt_landing_pads/fLP3/fLP3_Endo2_lb_expression_formatted.txt'

# infile <- '../../processed_data/endo_tss/alt_landing_pads/rLP6/rLP6_Endo2_lb_expression.txt'
# outfile <- '../../processed_data/endo_tss/alt_landing_pads/rLP6/rLP6_Endo2_lb_expression_formatted.txt'


data <- read.table(file = infile, header = T)

# parse name
data <- data %>% 
    mutate(name = gsub('>', '', orig_name),
           name = gsub('_rc', '', name)) %>% 
    separate(name, into = c('tss_name', 'tss_position', 'strand'), sep = ',', remove = F) %>% 
    mutate(tss_position = as.numeric(tss_position),
           start = ifelse(strand == '+', tss_position - 120, tss_position - 30),
           end = ifelse(strand == '+', tss_position + 30, tss_position + 120)) %>% 
    select(name:end, variant:num_barcodes_integrated, -orig_name)

data$category <- "tss"
data$category[grep("pos_control", data$name)] <- "pos_control"
data$category[grep("neg_control", data$name)] <- "neg_control"

write.table(data, file = outfile, quote = F, row.names = F)
