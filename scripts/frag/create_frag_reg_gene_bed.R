options(stringsAsFactors = F)

library(dplyr)
library(ggplot2)
require(cowplot)
library(sqldf)
library(matrixStats)
library(reshape2)

M9_num <- read.table("../../processed_data/frag/m9/Unshared_M9_peaks.bed") %>% nrow()

LB_num <- read.table("../../processed_data/frag/lb/Unshared_LB_peaks.bed") %>% nrow()

M9_total <- read.table("../../processed_data/frag/m9/U00096.2_M9_plus_minus_called_peaks_threshold1.1_merge40_min60.bed") %>% nrow()

shared <- M9_total - M9_num

# print(paste(M9_num, "unique to M9"))
# print(paste(LB_num, "unique to LB"))
# print(paste(shared, "shared"))

# Gene coordinates
gene_coords <- read.table("../../ref/U00096.2_genes_clean.bed", header = F) %>%
    select(gene = V4, Leftbound = V2, Rightbound = V3, strand = V6) %>%
    lapply(., gsub, pattern='NC_000913.2:', replacement='') %>%
    as.data.frame()

# Promoter peak with position of
promoter_max_LB <- read.table("../../processed_data/frag/lb/U00096.2_plus_minus_called_peaks_threshold1.1_merge40_min60.bed", header = F) %>%
    select(peak_name = V4, peak_max = V8, peak_score = V7, peak_strand =  V6) %>%
    mutate(downstream = ifelse(peak_strand == '+', peak_max + 500, peak_max-500), Condition = "LB") 

promoter_max_M9 <- read.table("../../processed_data/frag/m9/U00096.2_M9_plus_minus_called_peaks_threshold1.1_merge40_min60.bed", header = F) %>%
    select(peak_name = V4, peak_max = V8, peak_score = V7, peak_strand =  V6) %>%
    mutate(downstream = ifelse(peak_strand == '+', peak_max + 500, peak_max-500), Condition = "M9")


# match promoter to gene if contained within gene
intragenic_overlap_LB <- sqldf("SELECT a.*, b.*
                               FROM promoter_max_LB AS a 
                               LEFT JOIN gene_coords AS b
                               ON a.peak_max BETWEEN b.Leftbound AND b.Rightbound") %>% na.omit() %>%
    mutate(Orientation = ifelse(strand==peak_strand, "Sense", "Antisense"), 
           type = 'Intragenic')

intergenic_promoters_LB <- anti_join(promoter_max_LB, intragenic_overlap_LB, by = 'peak_name')

# Match promoter peak with nearest downstream gene
intergenic_overlap_LB <- sqldf("SELECT a.*, b.*
                               FROM intergenic_promoters_LB AS a 
                               JOIN gene_coords AS b
                               ON b.Leftbound BETWEEN a.peak_max AND a.downstream
                               OR b.Rightbound BETWEEN a.peak_max AND a.downstream
                               OR b.Leftbound BETWEEN a.downstream AND a.peak_max
                               OR b.Rightbound BETWEEN a.downstream AND a.peak_max") %>%
    na.omit() %>%
    group_by(peak_max) %>% 
    filter(if (peak_strand == '+') Leftbound == min(Leftbound) else Rightbound == max(Rightbound)) %>% 
    mutate(Orientation = ifelse(strand==peak_strand, "Sense", "Antisense"), 
           type =  'Intergenic') %>% ungroup()


# match promoter to gene if contained within gene
intragenic_overlap_M9 <- sqldf("SELECT a.*, b.*
                               FROM promoter_max_M9 AS a 
                               LEFT JOIN gene_coords AS b
                               ON a.peak_max BETWEEN b.Leftbound AND b.Rightbound") %>% na.omit() %>%
    mutate(Orientation = ifelse(strand==peak_strand, "Sense", "Antisense"), 
           type = 'Intragenic')


intergenic_promoters_M9 <- anti_join(promoter_max_M9, intragenic_overlap_M9, by = 'peak_name')

# Match non-intragenic promoter peak with nearest downstream gene
intergenic_overlap_M9 <- sqldf("SELECT a.*, b.*
                               FROM intergenic_promoters_M9 AS a 
                               JOIN gene_coords AS b
                               ON b.Leftbound BETWEEN a.peak_max AND a.downstream
                               OR b.Rightbound BETWEEN a.peak_max AND a.downstream
                               OR b.Leftbound BETWEEN a.downstream AND a.peak_max
                               OR b.Rightbound BETWEEN a.downstream AND a.peak_max") %>%
    na.omit() %>%
    group_by(peak_max) %>% 
    filter(if (peak_strand == '+') Leftbound == min(Leftbound) else Rightbound == max(Rightbound)) %>% 
    mutate(Orientation = ifelse(strand==peak_strand, "Sense", "Antisense"), 
           type =  'Intergenic') %>% ungroup()

# get sense reg genes
sense_reg_genes <- intergenic_overlap_M9 %>% filter(peak_strand == strand) %>% select(gene) %>% distinct()

# get antisense genes
temp <- intergenic_overlap_M9 %>%
    filter(peak_strand != strand) %>%
    select(gene)

antisense_reg_genes <- intragenic_overlap_M9 %>%
    filter(peak_strand != strand) %>%
    select(gene) %>%
    rbind(., temp) %>%
    distinct()


dual_reg_genes <- inner_join(antisense_reg_genes, sense_reg_genes, by = 'gene') %>% distinct()

# write out bed for plus dual reg
gene_coords %>% semi_join(dual_reg_genes, by = 'gene') %>%
    mutate(genome = 'U00096.2', x = '.') %>%
    select(genome, Leftbound, Rightbound, gene, x, strand) %>%
    filter(strand == '+') %>%
    write.table("../../processed_data/frag/m9/plus_dual_reg_genes.bed",
                col.names = FALSE,
                quote = F,
                row.names = FALSE,
                sep = '\t')

# Write out bed for minus dual reg
gene_coords %>% semi_join(dual_reg_genes, by = 'gene') %>%
    mutate(genome = 'U00096.2', x = '.') %>%
    select(genome, Leftbound, Rightbound, gene, x, strand) %>%
    filter(strand == '-') %>%
    write.table("../../processed_data/frag/m9/minus_dual_reg_genes.bed",
                col.names = FALSE,
                quote = F,
                row.names = FALSE,
                sep = '\t')

# Write out bed for plus sense reg
gene_coords %>% semi_join(sense_reg_genes, by = 'gene') %>% 
    anti_join(dual_reg_genes, by = 'gene') %>%
    mutate(genome = 'U00096.2', x = '.') %>%
    select(genome, Leftbound, Rightbound, gene, x, strand) %>%
    filter(strand == '+') %>%
    write.table("../../processed_data/frag/m9/plus_sense_reg_genes.bed",
                col.names = FALSE,
                quote = F,
                row.names = FALSE,
                sep = '\t')


# Write out bed for minus sense reg
gene_coords %>% semi_join(sense_reg_genes, by = 'gene') %>% 
    anti_join(dual_reg_genes, by = 'gene') %>%
    mutate(genome = 'U00096.2', x = '.') %>%
    select(genome, Leftbound, Rightbound, gene, x, strand) %>%
    filter(strand == '-') %>%
    write.table("../../processed_data/frag/m9/minus_sense_reg_genes.bed",
                col.names = FALSE,
                quote = F,
                row.names = FALSE,
                sep = '\t')

# Write out bed for plus antisense reg
gene_coords %>% semi_join(antisense_reg_genes, by = 'gene') %>% 
    anti_join(dual_reg_genes, by = 'gene') %>%
    mutate(genome = 'U00096.2', x = '.') %>%
    select(genome, Leftbound, Rightbound, gene, x, strand) %>%
    filter(strand == '+') %>%
    write.table("../../processed_data/frag/m9/plus_antisense_reg_genes.bed",
                col.names = FALSE,
                quote = F,
                row.names = FALSE,
                sep = '\t')


# Write out bed for minus antisense reg
gene_coords %>% semi_join(antisense_reg_genes, by = 'gene') %>% 
    anti_join(dual_reg_genes, by = 'gene') %>%
    mutate(genome = 'U00096.2', x = '.') %>%
    select(genome, Leftbound, Rightbound, gene, x, strand) %>%
    filter(strand == '-') %>%
    write.table("../../processed_data/frag/m9/minus_antisense_reg_genes.bed",
                col.names = FALSE,
                quote = F,
                row.names = FALSE,
                sep = '\t')