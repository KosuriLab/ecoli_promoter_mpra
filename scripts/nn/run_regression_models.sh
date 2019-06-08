# # define genome splits
# ./run_splits.sh regression

# # run neural network (works best on GPUs)
# ./run_nn.sh regression

# create k-mer counts
Rscript create_kmer_counts.R \
../../processed_data/combined/tss_scramble_peak_expression_model_format_floored_train_genome_split.txt \
../../processed_data/combined/tss_scramble_peak_expression_model_format_floored_train_genome_split_kmer_counts.txt

Rscript create_kmer_counts.R \
../../processed_data/combined/tss_scramble_peak_expression_model_format_floored_test_genome_split.txt \
../../processed_data/combined/tss_scramble_peak_expression_model_format_floored_test_genome_split_kmer_counts.txt