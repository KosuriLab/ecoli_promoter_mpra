OUTPUT=$1
DATE=`date +%Y%m%d`

echo "Running hyperparameter tuning for regression..."

python ../nn/dragonn_hyperparameter_tuning_regression.py \
peak_tile_expression_90train.txt \
peak_tile_expression_10test.txt \
150 4 5 100 0.2 100 $OUTPUT/$DATE_peak_tile_hyperparam_tuned > $OUTPUT/$DATE_tune_log.txt