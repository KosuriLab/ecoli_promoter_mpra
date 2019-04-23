DATA=$1
DATE=`date +%Y%m%d`

echo "Running hyperparameter tuning for classification..."

python ../nn/dragonn_hyperparameter_tuning.py \
$DATA/peak_tile_expression_90train_active.fasta \
$DATA/peak_tile_expression_90train_inactive.fasta \
$DATA/peak_tile_expression_10test_active.fasta \
$DATA/peak_tile_expression_10test_inactive.fasta \
150 4 5 100 0.2 100 $DATA/${DATE}_peak_tile_classifier_hyperparam_tuned > \
$DATA/${DATE}_peak_tile_classifier_hyperparam_tuned.log
