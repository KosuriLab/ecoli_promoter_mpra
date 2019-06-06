model_type=$1

if [ $model_type == 'regression' ]; then
	python define_genome_splits.py 0.75 4639675 \
	../../processed_data/combined/tss_scramble_peak_expression_model_format.txt \
	../../processed_data/combined/tss_scramble_peak_expression_model_format.txt --floor
fi

if [ $model_type == 'classification' ]; then
	python define_genome_splits.py 0.75 4639675 \
	../../processed_data/combined/tss_scramble_peak_expression_model_format.txt \
	../../processed_data/combined/tss_scramble_peak_expression_model_format.txt \
	--classification --neg_threshold 0.75 --pos_threshold 1.25
fi