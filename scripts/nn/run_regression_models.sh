DATE=`date +%Y%m%d`
# # define genome splits
# ./run_splits.sh regression

# # run neural network (works best on GPUs)
# ./run_nn.sh regression

# run MLP for 3 to 6-mers
python mlp_kmer.py \
../../processed_data/combined/tss_scramble_peak_expression_model_format_floored_train_genome_split.txt \
../../processed_data/combined/tss_scramble_peak_expression_model_format_floored_test_genome_split.txt \
../../processed_data/combined/${DATE}_mlp_tss_3to6mer_scramble_peak_predictions.txt 3 6 100

# # create 6-mer counts for random forest
# Rscript create_kmer_counts.R \
# ../../processed_data/combined/tss_scramble_peak_expression_model_format_floored_train_genome_split.txt \
# ../../processed_data/combined/tss_scramble_peak_expression_model_format_floored_train_genome_split_kmer_counts.txt

# Rscript create_kmer_counts.R \
# ../../processed_data/combined/tss_scramble_peak_expression_model_format_floored_test_genome_split.txt \
# ../../processed_data/combined/tss_scramble_peak_expression_model_format_floored_test_genome_split_kmer_counts.txt