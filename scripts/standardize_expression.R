# standardize expression across all datasets. Formatted expression files should be
# create first, e.g. ending with _expression_formatted.txt . Create a "standard curve" based
# on the wild-type sequences in both the TSS and scramble library. Use robust
# linear fit to mitigate effects of outliers on the standard linear fit


options(scipen = 10000)
options(stringsAsFactors = F)

library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)


tss <- read.table('../processed_data/endo_tss/lb/rLP5_Endo2_lb_expression_formatted.txt', 
                  header = T)

scramble <- read.table('../processed_data/endo_scramble/endo_scramble_expression_formatted.txt', 
                       header = T)

peak_tile <- read.table('../processed_data/peak_tile/peak_tile_expression_formatted.txt',
                        header = T)

flp3 <- read.table('../processed_data/endo_tss/alt_landing_pads/fLP3/fLP3_Endo2_lb_expression_formatted.txt',
                   header = T)

rlp6 <- read.table('../processed_data/endo_tss/alt_landing_pads/rLP6/rLP6_Endo2_lb_expression_formatted.txt',
                   header = T)

# format scramble negative control category correctly
scramble$category[grep("_scrambled", scramble$name)] <- "scramble"

# wt_scramble <- scramble %>% 
#     filter(category == 'unscrambled') %>% 
#     left_join(select(tss, tss_name, RNA_exp_sum_ave, expn_med), by = 'tss_name',
#               suffix = c('_scramble', '_tss'))
# 
# wt_robust_fit <- MASS::rlm(expn_med_tss ~ expn_med_scramble, wt_scramble)
# 
# # fit expression using robust fit, rename column so it works with fit
# scramble$expn_med_fitted <- predict(wt_robust_fit, 
#                                     select(scramble, expn_med_scramble = expn_med))
# peak_tile$expn_med_fitted <- predict(wt_robust_fit,
#                                      select(peak_tile, expn_med_scramble = expn_med))
# flp3$expn_med_fitted <- predict(wt_robust_fit,
#                                      select(flp3, expn_med_scramble = expn_med))
# rlp6$expn_med_fitted <- predict(wt_robust_fit,
#                                      select(rlp6, expn_med_scramble = expn_med))

standardize <- function(df1, df2) {
    # subset positive controls
    pos_controls <- inner_join(filter(df1, category == 'pos_control') %>% 
                                   select(name, expn_med1 = expn_med),
                               filter(df2, category == 'pos_control') %>% 
                                   select(name, expn_med2 = expn_med),
                               by = 'name')
    control_fit <- MASS::rlm(expn_med1 ~ expn_med2, pos_controls)
    predicted <- predict(control_fit, select(df2, expn_med2 = expn_med))
    return(predicted)
}

scramble$expn_med_fitted <- standardize(tss, scramble)
peak_tile$expn_med_fitted <- standardize(tss, peak_tile)
flp3$expn_med_fitted <- standardize(tss, flp3)
rlp6$expn_med_fitted <- standardize(tss, rlp6)

# we do not need to scale TSS because it's the standard, but create a "fitted"
# column for consistency
tss$expn_med_fitted <- tss$expn_med

set_std_threshold <- function(df) {
    neg <- filter(df, category == 'neg_control')
    neg_sd <- sd(neg$expn_med_fitted)
    neg_median <- median(neg$expn_med_fitted)
    # threshold <- neg_median + (2 * neg_sd)
    # scale = 1 / threshold
    # print(c(threshold, scale))
    # return(df$expn_med_fitted * scale)
    return(df$expn_med_fitted - neg_median)
}

set_std_threshold_peak <- function(df) {
    neg <- filter(df, category == 'neg_control')
    neg_mad <- mad(neg$expn_med_fitted)
    neg_median <- median(neg$expn_med_fitted)
    # threshold <- neg_median + (3 * neg_mad)
    # scale = 1 / threshold
    # print(c(threshold, scale))
    # return(df$expn_med_fitted * scale)
    return(df$expn_med_fitted - neg_median)
}

tss$expn_med_fitted_scaled <- set_std_threshold(tss)
scramble$expn_med_fitted_scaled <- set_std_threshold(scramble)
peak_tile$expn_med_fitted_scaled <- set_std_threshold_peak(peak_tile)
flp3$expn_med_fitted_scaled <- set_std_threshold(flp3)
rlp6$expn_med_fitted_scaled <- set_std_threshold(rlp6)

# now that expression is standardized, calculate change in expression for scramble
# calculate change in expression: unscramble - scramble
scramble <- scramble %>% 
    group_by(tss_name) %>% 
    mutate(unscrambled_exp = ifelse(any(category == 'unscrambled'),
                                    RNA_exp_sum_ave[category == 'unscrambled'],
                                    NA),
           relative_exp = RNA_exp_sum_ave / unscrambled_exp)

# create active and inactive columns for TSS and alternate landing pads
tss <- tss %>% 
    mutate(active = ifelse(expn_med_fitted_scaled >= 1, 'active', 'inactive'))
peak_tile <- peak_tile %>% 
    mutate(active = ifelse(expn_med_fitted_scaled >= 1, 'active', 'inactive'))
flp3 <- flp3 %>% 
    mutate(active = ifelse(expn_med_fitted_scaled >= 1, 'active', 'inactive'))
rlp6 <- rlp6 %>% 
    mutate(active = ifelse(expn_med_fitted_scaled >= 1, 'active', 'inactive'))

write.table(tss, file = '../processed_data/endo_tss/lb/rLP5_Endo2_lb_expression_formatted_std.txt',
            row.names = F, quote = F, sep = '\t')
write.table(scramble, file = '../processed_data/endo_scramble/endo_scramble_expression_formatted_std.txt',
            row.names = F, quote = F, sep = '\t')
write.table(peak_tile, file = '../processed_data/peak_tile/peak_tile_expression_formatted_std.txt',
            row.names = F, quote = F, sep = '\t')
write.table(flp3, file = '../processed_data/endo_tss/alt_landing_pads/fLP3/fLP3_Endo2_lb_expression_formatted_std.txt',
            row.names = F, quote = F, sep = '\t')
write.table(rlp6, file = '../processed_data/endo_tss/alt_landing_pads/rLP6/rLP6_Endo2_lb_expression_formatted_std.txt',
            row.names = F, quote = F, sep = '\t')

# combine datasets for modeling
combined <- bind_rows(tss, 
                      select(scramble, variant, expn_med_fitted_scaled,
                             start=var_left, end=var_right, name),
                      peak_tile) %>% 
    select(variant, expn_med_fitted_scaled, start, end, name)

write.table(combined, '../processed_data/combined/tss_scramble_peak_expression_model_format.txt',
            row.names = F, col.names = F, quote = F, sep = '\t')
write.table(select(combined, variant, expn_med_fitted_scaled), 
            '../processed_data/combined/tss_scramble_peak_expression_model_format_values_only.txt',
            row.names = F, col.names = F, quote = F, sep = '\t')
# just TSS for gkmSVM
tss %>% 
    select(variant, expn_med_fitted_scaled, start, end, name) %>% 
    write.table('../processed_data/endo_tss/lb/model_files/tss_expression_model_format.txt',
                row.names = F, col.names = F, quote = F, sep = '\t')


# graph positive controls in scramble and TSS
pos_controls <- inner_join(filter(tss, category == 'pos_control') %>% 
                               select(name, expn_med1 = expn_med_fitted_scaled),
                           filter(scramble, category == 'pos_control') %>% 
                               select(name, expn_med2 = expn_med_fitted_scaled),
                           by = 'name')

controls <- inner_join(filter(tss, grepl("control", category)) %>% 
                           select(name, category, expn_med1 = expn_med_fitted),
                           filter(scramble, grepl("control", category)) %>% 
                               select(name, category, expn_med2 = expn_med_fitted),
                           by = 'name')

ggplot(controls, aes(expn_med1, expn_med2)) + 
    geom_point(aes(color = category.x)) +
    scale_x_log10() + scale_y_log10()

# graph unscrambled
wildtype <- inner_join(filter(tss, category == 'tss') %>% 
                           select(name, tss_name, expn_med1 = expn_med),
                       filter(scramble, category == 'unscrambled') %>% 
                           select(name, tss_name, expn_med2 = expn_med),
                       by = 'tss_name')


r2 <- summary(lm(expn_med1 ~ expn_med2, wildtype))$r.squared
ggplot(wildtype, aes(expn_med1, expn_med2)) + 
    geom_point() + 
    annotation_logticks(sides='bl') + scale_x_log10() + scale_y_log10() +
    geom_smooth(method = 'lm') +
    labs(x = 'TSS expression', y = 'scramble expression',
         title = 'Wild-type sequences') +
    annotate('text', x=50, y=0.1, parse=T, label=paste('R^2==', signif(r2, 3)))

