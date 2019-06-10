# # define genome splits
# ./run_splits.sh classification

# # run neural network (works best on GPUs)
# ./run_nn.sh classification

# create k-mer counts
Rscript create_kmer_counts.R \
../../processed_data/combined/tss_scramble_peak_expression_model_format_train_genome_split_classification.txt \
../../processed_data/combined/tss_scramble_peak_expression_model_format_train_genome_split_classification_kmer_counts.txt

Rscript create_kmer_counts.R \
../../processed_data/combined/tss_scramble_peak_expression_model_format_test_genome_split_classification.txt \
../../processed_data/combined/tss_scramble_peak_expression_model_format_test_genome_split_classication_kmer_counts.txt

