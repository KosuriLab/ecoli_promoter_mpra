library(dplyr)
library(tidyr)

options(stringsAsFactors = F)

args = commandArgs(trailingOnly=TRUE)

# infile <- args[1]
# outfile <- args[2]
infile <- '../../processed_data/endo_tss/lb/rLP5_Endo2_lb_expression.txt'
outfile <- '../../processed_data/endo_tss/lb/rLP5_Endo2_lb_expression_formatted.txt'

data <- read.table(file = infile, header = T)

# parse name
data <- data %>% 
    mutate(name = gsub('>', '', orig_name),
           name = gsub('_rc', '', name)) %>% 
    separate(name, into = c('source', 'tss_pos', 'strand'), sep = ',', remove = F) %>% 
    select(name:strand, variant:num_barcodes_integrated, -orig_name)

data$category <- "tss"
data$category[grep("pos_control", data$name)] <- "pos_control"
data$category[grep("neg_control", data$name)] <- "neg_control"

write.table(data, file = outfile, quote = F, row.names = F)
