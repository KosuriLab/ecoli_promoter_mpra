options(stringsAsFactors = F)
library(dplyr)
library(tidyr)


data <- read.table('../../processed_data/peak_tile/peak_tile_expression_formatted_std.txt',
                   header = T)

# split into pre-define train and test based on peak
data <- data %>% 
    unite(peak_name, peak_start, peak_end, sep = 'to', remove = F)
peak_names <- unique(data$peak_name)

set.seed(123)
# get 90% of peaks
train_size <- floor(0.90 * length(peak_names))
train_index <- sample(seq_len(length(peak_names)), size = train_size)
train_peak_names <- peak_names[train_index]
test_peak_names <- peak_names[-train_index]
data_train <- filter(data, peak_name %in% train_peak_names)
data_test <- filter(data, peak_name %in% test_peak_names)

write.table(select(data_train, variant, expn_med_fitted_scaled),
            file = '../../processed_data/peak_tile/peak_tile_expression_90train.txt',
            sep  = '\t', col.names = F, row.names = F, quote = F)

write.table(select(data_test, variant, expn_med_fitted_scaled),
            file = '../../processed_data/peak_tile/peak_tile_expression_10test.txt',
            sep  = '\t', col.names = F, row.names = F, quote = F)
