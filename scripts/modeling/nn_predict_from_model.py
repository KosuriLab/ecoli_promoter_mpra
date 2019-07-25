# from __future__ import absolute_import, division, print_function
import numpy as np, random
np.random.seed(1)
random.seed(1)
from keras_regression import SequenceDNN
from keras.models import load_model
from hyperparameter_search_regression import HyperparameterSearcher, RandomSearch
from sklearn.preprocessing import LabelEncoder, OneHotEncoder
try:
    from sklearn.model__selection import train_test_split  # sklearn >= 0.18
except ImportError:
    from sklearn.cross_validation import train_test_split  # sklearn < 0.18
import sys
import argparse


num_epochs = 100


def one_hot_encode(sequences):
	# horizontal one-hot encoding
    sequence_length = len(sequences[0])
    integer_type = np.int8 if sys.version_info[
        0] == 2 else np.int32  # depends on Python version
    integer_array = LabelEncoder().fit(np.array(('ACGTN',)).view(integer_type)).transform(
        sequences.view(integer_type)).reshape(len(sequences), sequence_length)
    one_hot_encoding = OneHotEncoder(
        sparse=False, n_values=5, dtype=integer_type).fit_transform(integer_array)

    return one_hot_encoding.reshape(
        len(sequences), 1, sequence_length, 5).swapaxes(2, 3)[:, :, [0, 1, 2, 4], :]


def pad_sequence(seq, max_length):
    if len(seq) > max_length:
        diff = len(seq) - max_length
        # diff%2 returns 1 if odd
        trim_length = int(diff / 2)
        seq = seq[trim_length : -(trim_length + diff%2)]
    else:
        seq = seq.center(max_length, 'N')

    return seq


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('sequences', help='''tab-separated, two columns. 
        First is sequence, second is continuous value for regression or 0/1 for
        classification''')
    parser.add_argument('seq_length', type=int, help='length of input sequences')
    parser.add_argument('arch_file', help='.json of model architecture')
    parser.add_argument('weights_file', help='h5 of model weights')
    parser.add_argument('prediction_file', help='output name of predictions')

    args = parser.parse_args()

    sequences = args.sequences
    seq_length = args.seq_length
    arch_file = args.arch_file
    weights_file = args.weights_file

    # read in sequences and labels
    print("loading sequence data...")
    seqs = [line.split('\t')[0] for line in open(sequences)]
    padded_seqs = [pad_sequence(x, seq_length) for x in seqs]
    X = one_hot_encode(np.array(padded_seqs))
    y = np.array([float(line.strip().split('\t')[1]) for line in open(sequences)])

    model = SequenceDNN.load(arch_file, weights_file)
    predictions = np.squeeze(model.predict(X))
    corr = np.corrcoef(predictions, y)[0,1]
    print("Test set correlation:", corr)
    # write predictions
    with open(args.prediction_file, 'w') as outfile:
        outfile.write('observed\tpredicted\n')
        for i in range(len(y)):
            outfile.write(str(y[i]) + '\t' + str(predictions[i]) + '\n')



    