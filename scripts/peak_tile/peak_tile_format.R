options(scipen = 10000)
options(stringsAsFactors = F)

library(dplyr)
library(tidyr)

data <- read.table('../../processed_data/peak_tile/peak_tile_expression.txt', header = T) %>% 
    select(variant, name = orig_name, RNA_exp_sum_1_1:num_barcodes_integrated)

data  <- data  %>% 
    mutate(name = gsub('_flipped', '', name),
           name = gsub('_rc', '', name))

data <- data%>% 
    mutate(category = case_when(grepl('random', .$name) ~ 'random',
                                grepl('neg', .$name) ~ 'neg_control',
                                grepl('pos_control', .$name) ~ 'pos_control',
                                TRUE ~ 'tile'))

data <- data %>% 
    separate(name, into = c('peak_start', 'peak_end', 'strand', 'tile_loc'), 
             sep = '_', remove = F, convert = T) %>% 
    mutate(peak_start = as.numeric(peak_start),
           peak_end = as.numeric(peak_end),
           peak_length = peak_end - peak_start)

data <- data %>% 
    mutate(tile_loc = gsub('pos', '', tile_loc)) %>% 
    separate(tile_loc, into = c('tile_start', 'tile_end'),
             sep = '-', convert = T) %>% 
    mutate(tile_start = as.numeric(tile_start),
           tile_start_relative = tile_start/peak_length)

data <- data %>% 
    mutate(start = peak_start + tile_start,
           end = peak_end + tile_end)

write.table(data, file = '../../processed_data/peak_tile/peak_tile_expression_formatted.txt',
            quote = F, row.names = F, sep = '\t')
