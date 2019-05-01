DATA=$1
DATE=`date +%Y%m%d`

echo "Running hyperparameter tuning for regression..."

# touch ${DATA}/${DATE}_peak_tile_hyperparam_tuned.log
python ../nn/dragonn_hyperparameter_tuning_regression.py \
$DATA/peak_tile_expression_90train.txt \
$DATA/peak_tile_expression_10test.txt \
150 4 5 100 0.2 100 $DATA/model_files/${DATE}_peak_tile_hyperparam_tuned > \
$DATA/model_files/${DATE}_peak_tile_hyperparam_tuned.log
