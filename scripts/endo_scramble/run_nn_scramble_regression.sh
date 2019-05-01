TRAIN=$1
TEST=$2
OUTPUT=$3
DATE=`date +%Y%m%d`

echo "Running hyperparameter tuning for regression..."

# touch ${DATA}/${DATE}_peak_tile_hyperparam_tuned.log
python ../nn/dragonn_hyperparameter_tuning_regression.py \
$TRAIN $TEST \
150 4 5 100 0.2 100 $OUTPUT/${DATE}_scramble_hyperparam_tuned > $OUTPUT/${DATE}_scramble_hyperparam_tuned.log
