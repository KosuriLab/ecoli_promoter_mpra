model_type=$1
DATE=`date +%Y%m%d`

if [ $model_type == 'regression' ]; then
	python dragonn_hyperparameter_tuning_regression.py \
	../../processed_data/combined/tss_scramble_peak_expression_model_format_floored_train_genome_split.txt \
	../../processed_data/combined/tss_scramble_peak_expression_model_format_floored_test_genome_split.txt \
	150 4 5 100 0.2 100 \
	../../processed_data/combined/${DATE}_tss_scramble_peak_regression_hyperparam_tuned > \
	../../processed_data/combined/${DATE}_tss_scramble_peak_regression_hyperparam_tuned.log
fi

if [ $model_type == 'classification' ]; then
	python dragonn_hyperparameter_tuning.py \
	../../processed_data/combined/tss_scramble_peak_expression_model_format_train_genome_split_classification.txt \
	../../processed_data/combined/tss_scramble_peak_expression_model_format_test_genome_split_classification.txt \
	150 4 5 100 0.2 100 \
	../../processed_data/combined/${DATE}_tss_scramble_peak_classification_hyperparam_tuned > \
	../../processed_data/combined/${DATE}_tss_scramble_peak_classification_hyperparam_tuned.log
fi
