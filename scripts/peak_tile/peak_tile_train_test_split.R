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


writeFasta <- function(data, filename){
    fastaLines = c()
    for (rowNum in 1:nrow(data)){
        fastaLines = c(fastaLines, as.character(paste(">", data[rowNum,"name"], sep = "")))
        fastaLines = c(fastaLines,as.character(data[rowNum,"variant"]))
    }
    fileConn<-file(filename)
    writeLines(fastaLines, fileConn)
    close(fileConn)
}

# write files for classification, split by active/inactive
writeFasta(filter(data_train, active == 'active') %>% select(name, variant),
           filename = '../../processed_data/peak_tile/peak_tile_expression_90train_active.fasta')

writeFasta(filter(data_train, active == 'inactive') %>% select(name, variant),
           filename = '../../processed_data/peak_tile/peak_tile_expression_90train_inactive.fasta')

writeFasta(filter(data_test, active == 'active') %>% select(name, variant),
           filename = '../../processed_data/peak_tile/peak_tile_expression_10test_active.fasta')

writeFasta(filter(data_test, active == 'inactive') %>% select(name, variant),
           filename = '../../processed_data/peak_tile/peak_tile_expression_10test_inactive.fasta')
