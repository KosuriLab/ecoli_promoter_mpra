# from __future__ import absolute_import, division, print_function
import numpy as np, random
np.random.seed(1)
random.seed(1)
from keras_regression import SequenceDNN
from hyperparameter_search_regression import HyperparameterSearcher, RandomSearch
from sklearn.preprocessing import LabelEncoder, OneHotEncoder
try:
    from sklearn.model_selection import train_test_split  # sklearn >= 0.18
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


def one_hot_encode_2d(sequences):
	# horizontal one-hot encoding
    sequence_length = len(sequences[0])
    integer_type = np.int8 if sys.version_info[
        0] == 2 else np.int32  # depends on Python version
    integer_array = LabelEncoder().fit(np.array(('ACGTN',)).view(integer_type)).transform(
        sequences.view(integer_type)).reshape(len(sequences), sequence_length)
    one_hot_encoding = OneHotEncoder(
        sparse=False, n_values=5, dtype=integer_type).fit_transform(integer_array)
    # dimensions are n-samples, n-features. The one hot encoded vector is kept as a single
    # vector instead of split into 1x4 matrix. n-features = 4 * sequence_length
    return one_hot_encoding


def pad_sequence(seq, max_length):
	if len(seq) > max_length:
		diff = len(seq) - max_length
		# diff%2 returns 1 if odd
		trim_length = int(diff / 2)
		seq = seq[trim_length : -(trim_length + diff%2)]
	else:
		seq = seq.center(max_length, 'N')

	return seq


def process_seqs(filename, seq_length, encode_type='1d'):

	seqs = [line.split('\t')[0] for line in open(filename)]
	padded_seqs = [pad_sequence(x, seq_length) for x in seqs]
	if encode_type == '1d':
		X = one_hot_encode(np.array(padded_seqs))
	elif encode_type == '2d':
		# only keep 150bp sequences
		X = one_hot_encode_2d(np.array([x for x in padded_seqs]))
	else:
		raise ValueException('Non-valid encoding type')

	y = np.array([float(line.strip().split('\t')[1]) for line in open(filename)])
	return X, y


if __name__ == '__main__':

	parser = argparse.ArgumentParser()
	parser.add_argument('data', help='''Training set, one column sequence, 
		one column expression. Tab-separated''')
	parser.add_argument('seq_length', type=int, help='length of input sequences')
	args = parser.parse_args()

	# load in pre-defined splits
	seq_length = args.seq_length
	print("loading training set...")
	X, y = process_seqs(args.train, seq_length, '2d')

	model = keras_regression.SVR()
	print("Training SVR...")
	start_time = time.time()
	model.train(X, y)
	elapsed_time = time.time() - start_time
	print(elapsed_time)
	print "Score:"
	model.score(X, y)
	










