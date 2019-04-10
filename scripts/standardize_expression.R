# standardize expression across all datasets. Formatted expression files should be
# create first, e.g. ending with _expression_formatted.txt . Create a "standard curve" based
# on the wild-type sequences in both the TSS and scramble library. Use robust
# linear fit to mitigate effects of outliers on the standard linear fit


options(scipen = 10000)
options(stringsAsFactors = F)

library(dplyr)
library(tidyr)


tss <- read.table('../processed_data/endo_tss/lb/rLP5_Endo2_lb_expression_formatted.txt', 
                  header = T)

scramble <- read.table('../processed_data/endo_scramble/endo_scramble_expression_formatted.txt', 
                       header = T)

peak_tile <- read.table('../processed_data/peak_tile/peak_tile_expression_formatted.txt',
                        header = T)

wt_scramble <- scramble %>% 
    filter(category == 'unscrambled') %>% 
    left_join(select(tss, tss_name, RNA_exp_sum_ave, expn_med), by = 'tss_name',
              suffix = c('_scramble', '_tss'))

wt_robust_fit <- MASS::rlm(expn_med_tss ~ expn_med_scramble, wt_scramble)

# fit expression using robust fit, rename column so it works with fit
scramble$expn_med_fitted <- predict(wt_robust_fit, 
                                    select(scramble, expn_med_scramble = expn_med))
peak_tile$expn_med_fitted <- predict(wt_robust_fit,
                                     select(peak_tile, expn_med_scramble = expn_med))

# we do not need to scale TSS because it's the standard, but create a "fitted"
# column for consistency
tss$expn_med_fitted <- tss$expn_med

set_std_threshold <- function(df) {
    neg <- filter(df, category == 'neg_control')
    neg_sd <- sd(neg$expn_med_fitted)
    neg_median <- median(neg$expn_med_fitted)
    threshold <- neg_median + (2 * neg_sd)
    scale = 1 / threshold
    return(df$expn_med_fitted * scale)
}

tss$expn_med_fitted_scaled <- set_std_threshold(tss)
scramble$expn_med_fitted_scaled <- set_std_threshold(scramble)
peak_tile$expn_med_fitted_scaled <- set_std_threshold(peak_tile)

# now that expression is standardized, calculate change in expression for scramble
# calculate change in expression: unscramble - scramble
scramble <- scramble %>% 
    group_by(tss_name) %>% 
    mutate(unscrambled_exp = ifelse(any(category == 'unscrambled'),
                                    RNA_exp_sum_ave[category == 'unscrambled'],
                                    NA),
           relative_exp = RNA_exp_sum_ave / unscrambled_exp)

write.table(tss, file = '../processed_data/endo_tss/lb/rLP5_Endo2_lb_expression_formatted_std.txt',
            row.names = F, quote = F, sep = '\t')
write.table(scramble, file = '../processed_data/endo_scramble/endo_scramble_expression_formatted_std.txt',
            row.names = F, quote = F, sep = '\t')
write.table(scramble, file = '../processed_data/peak_tile/peak_tile_expression_formatted_std.txt',
            row.names = F, quote = F, sep = '\t')

