library(dplyr)
library(tidyr)
library(stringr)
options(stringsAsFactors = F)

data <- read.table('../../ref/20180508_lb_peak_tile_lib.txt', col.names = c('name', 'sequence'))

data <- data %>% 
    separate(name, into = c('peak_start', 'peak_end', 'strand', 'tile_loc'), 
             sep = '_', remove = F, convert = T) %>% 
    mutate(peak_start = as.numeric(peak_start),
           peak_end = as.numeric(peak_end),
           peak_length = peak_end - peak_start,
           variant = substr(sequence, 25, 174))

short_peaks <- filter(data, peak_length < 150)

# remove leading stuffer of T's
short_seqs <- str_split(short_peaks$variant, '^T+')
# select correct element from sublist, if short peak is long enough then it will not have a stuffer
# and there is only one element in sublist
short_seqs <- unlist(short_seqs)
short_seqs <- short_seqs[short_seqs != ""]
short_peaks$trimmed_seq <- short_seqs

# write file
write.table(select(short_peaks, name, trimmed_seq),
            '../../ref/20180508_lb_peak_tile_lib_short_peaks_trimmed.txt',
            sep = '\t', row.names = F, quote = F, col.names = F)
