# This script combines TSS information from three different sources: 
# 14,865 putative TSS from Thomason and Storz, 2,122 TSS from Wanner and 8,480 
# TSS from RegulonDB. This list is then collapse into a unique set of TSS which
# cannot be closer than 20bp to each other.
#
# Author: Kimberly Insigne, kiminsigne@gmail.com

library(dplyr)
options(stringsAsFactors=F)

tss_map <- read.table('tss_map_storz.txt', header=T, comment.char='#', sep='\t')

# Promoters can be different classes depending on the condition, so condition
# information is included for this study but is not relevant in the other two
# TSS sets

# In this study, promoters can be one of multiple classes: primary, secondary
# internal or antisense. Represent this information as 0 if it not this class,
# 1 otherwise. Combine all this information into one field like 0,1,0,0 using 
# the paste command
tss1 <- data.frame(tss_map$Pos, tss_map$Condition, tss_map$Strand, 
                 paste(tss_map$Primary, tss_map$Secondary, tss_map$Internal, 
                       tss_map$Antisense, sep=','))
colnames(tss1) <- c('position', 'condition', 'strand', 'type')
# add column for source information
tss1 <- cbind(tss1, source=rep('storz', nrow(tss1)))

# read in other tss file
tss2 <- read.table('tss_only_wanner.txt', header=F, sep='\t')
colnames(tss2) <- c('position', 'strand', 'type')
# Put in blank condition column
tss2 <- cbind(tss2, condition=rep(NA, nrow(tss2)))
# add source information
tss2 <- cbind(tss2, source=rep('wanner', nrow(tss2)))

# read in regulondb promoters
regulondb <- read.table('regulondb_promoters.txt', comment.char='#', 
                        header=F, sep='\t')
colnames(regulondb) <- c('id', 'name', 'strand', 'position', 'sigma_factor', 
                        'sequence', 'evidence')

# change forward and reverse strands to '+' and '-'
regulondb <- regulondb %>%
    mutate(strand = ifelse(strand == 'forward', '+', '-'))

# extract info to add to final_tss, add blank columns for condition and type
regulon_tss <- data.frame(position = regulondb$position, strand = regulondb$strand, 
            condition=rep(NA, nrow(regulondb)), type=rep(NA, nrow(regulondb)), 
            source=rep('regulondb', nrow(regulondb)))

# get rid of the 0's, don't care about these
regulon_tss <- filter(regulon_tss, position != 0)

# bind all TSSs together, sort by position
final_tss <- bind_rows(tss1, tss2, regulon_tss) %>% 
    arrange(position)

# how many TSSs shared by all three sources?
tss_in_all_three <- final_tss %>% 
    select(-condition, -type) %>% 
    distinct() %>% 
    group_by(position) %>% 
    summarise(num_sources = n()) %>% 
    filter(num_sources == 3)

# grab list of unique TSS and strand information
tss_only <- select(final_tss, position, strand) %>% distinct()

# for each unique TSS, create new column with list of all sources
final_tss <- final_tss %>%
    group_by(strand, position) %>%
    mutate(sources = paste(unique(source), collapse='_')) %>%
    ungroup()

# trim down unnecessary columns, grab unique
final_tss_unique <- final_tss %>%
    select(position, strand, -source, source=sources) %>%
    distinct()

# Takes a list of TSSs and builds up a list of TSSs, starting from smallest, so that no TSSs are within
# n bp of each other
distinct_tss <- function(tss_df, n){
    if(length(unique(tss_df$strand)) != 1){
        stop('More than one strand. Please separate so that only one strand is present.')
    }
    tss <- unlist(tss_df$position)
    # membership vector used for subsetting
    include <- rep(FALSE, length(tss))
    # initialize
    prev_tss <- tss[1]
    include[1] <- TRUE
    for(i in 2:length(tss)){
        curr_tss <- tss[i]
        if(curr_tss - prev_tss >= n){
            include[i] <- TRUE
            prev_tss <- curr_tss
        }
    }
    return(tss_df[include,])
}

# trim down so no TSSs are closer than 20bp, strand-dependent manner
forward_distinct <- distinct_tss(filter(final_tss_unique, strand == '+'), 20)
reverse_distinct <- distinct_tss(filter(final_tss_unique, strand == '-'), 20)

tss_trimmed <- bind_rows(forward_distinct, reverse_distinct) %>%
    arrange(position)

# write to output
write.table(tss_trimmed, file='final_tss_list.txt', quote=F, 
            row.names=F, sep='\t')


# make bar graph showing amount of overlap between different sources
# count how many are in all three sources
# count number of sources for each TSS
df <- select(final_tss, -condition) %>%
    distinct(position, strand, source)

# df <- df %>%
#     group_by(position, strand) %>%
#     summarize(num_sources = n())

# each soure has own column
df <- df %>%
    mutate(storz = ifelse(source == 'storz', 1, 0),
           wanner = ifelse(source == 'wanner', 1, 0),
           regulon = ifelse(source == 'regulondb', 1, 0)) %>%
    select(-source)

df <- df %>%
    group_by(position, strand) %>%
    summarize(storz = sum(storz), wanner = sum(wanner), regulon = sum(regulon))

# difference between consecutive TSSs
pos <- filter(df, strand == '+')
neg <- filter(df, strand == '-')
pos$diff <- c(Inf, diff(pos$position))
neg$diff <- c(Inf, diff(neg$position))
df <- bind_rows(pos, neg) %>%
    arrange(position)

df <- df %>%
    mutate(num_sources = storz + wanner + regulon)
# copy
unique_df <- df

# for(i in 2:nrow(unique_df)){
#     if(unique_df$diff[i] <= 3){
#         unique_df$num_sources[i] <- max(unique_df$num_sources[i-1], unique_df$num_sources[i])
#         unique_df$num_sources[i-1] <- max(unique_df$num_sources[i-1], unique_df$num_sources[i])
#     }
# }

for(i in 2:nrow(unique_df)){
    if(unique_df$diff[i] <= 3){
        total_source <- sum(unique_df$storz[(i-1):i]) + sum(unique_df$wanner[(i-1):i]) + sum(unique_df$regulon[(i-1):i])
        unique_df$num_sources[i] <- total_source
        unique_df$num_sources[i-1] <- total_source
    }
}

# get rid of positions that had an adjacent TSS <= 3
unique_df <- filter(unique_df, diff > 3)

unique_df <- unique_df %>%
    mutate(num_sources = ifelse(num_sources > 3, 3, num_sources))


labels <- unique_df %>%
    group_by(num_sources) %>%
    summarize(percent = round(sum(num_sources)/nrow(unique_df) * 100, 2))
labels$lab <- as.character(labels$percent)

ggplot(unique_df, aes(num_sources)) + 
    geom_bar(aes(y = (..count..)/sum(..count..) * 100), fill='navyblue') +
    geom_text(data = labels ,aes(x = num_sources, y = percent, label = lab), 
              vjust = 0, size = 6) +
    labs(x='Number of sources', 
         y='percent of TSSs',
         title='Overlap in TSSs reported from\ntwo genome-wide studies and RegulonDB')

