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
    integer_array = LabelEncoder().fit(np.array(('ACGT',)).view(integer_type)).transform(
        sequences.view(integer_type)).reshape(len(sequences), sequence_length)
    one_hot_encoding = OneHotEncoder(
        sparse=False, n_values=4, dtype=integer_type).fit_transform(integer_array)
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


def process_seqs(filename, seq_length):

	seqs = [line.split('\t')[0] for line in open(filename)]
	padded_seqs = [pad_sequence(x, seq_length) for x in seqs]
	X = one_hot_encode(np.array(padded_seqs))
	y = np.array([float(line.strip().split('\t')[1]) for line in open(filename)])
	return X, y


if __name__ == '__main__':

	parser = argparse.ArgumentParser()
	parser.add_argument('train', help='''pre-defined training set, 
		one column sequence, one column expression. Tab-separated''')
	parser.add_argument('test', help='pre-defined test set')
	parser.add_argument('seq_length', type=int, help='length of input sequences')
	parser.add_argument('num_layers', type=int, help='number of convolutional layers')
	parser.add_argument('min_filter', type=int, help='minimum number of filters')
	parser.add_argument('max_filter', type=int, help='maximum number of filters')
	parser.add_argument('validation_fraction', type=float)
	parser.add_argument('num_trials', type=int, 
		help='number of hyperparameter trials')
	parser.add_argument('prefix', help='output prefix for saved model files')
	
	args = parser.parse_args()

	# load in pre-defined splits
	seq_length = args.seq_length
	print("loading training set...")
	X_train, y_train = process_seqs(args.train, seq_length)
	print("loading test set...")
	X_test, y_test = process_seqs(args.test, seq_length)

	
	num_layers = args.num_layers
	min_filter = args.min_filter
	max_filter = args.max_filter
	validation_fraction = args.validation_fraction
	num_hyperparameter_trials = args.num_trials
	prefix = args.prefix


	# split training into validation set
	X_train, X_valid, y_train, y_valid = train_test_split(X_train, y_train, 
		test_size=validation_fraction)


	print('Starting hyperparameter search...')
	min_layer = 1
	max_layer = 4
	min_conv_width = 15
	max_conv_width = 20
	min_dropout = 0.1
	max_dropout = 0.9

	fixed_hyperparameters = {'seq_length': seq_length, 'num_epochs': num_epochs}
	grid = {'num_filters': ((min_filter, max_filter),), 'pool_width': (5, 40),
	        'conv_width': ((min_conv_width, max_conv_width),), 
	        'dropout': (min_dropout, max_dropout)}

	# number of convolutional layers        
	print("Number of convolutional layers: ", num_layers)
	filters = tuple([(min_filter, max_filter)] * num_layers)
	conv_widths = tuple([(min_conv_width, max_conv_width)] * num_layers)
	grid.update({'num_filters': filters, 'conv_width': conv_widths})

	# Backend is RandomSearch; if using Python 2, can also specify MOESearch
	# (requires separate installation)
	searcher = HyperparameterSearcher(SequenceDNN, fixed_hyperparameters, grid, 
		X_train, y_train, validation_data=(X_valid, y_valid), backend=RandomSearch)
	searcher.search(num_hyperparameter_trials)
	print('Best hyperparameters: {}'.format(searcher.best_hyperparameters))
	model = searcher.best_model
	
	# print test results
	print('Test results: {}'.format(model.score(X_test, y_test)))
	
	# save model
	model.save(prefix)
	
	# print predictions
	predictions = model.predict(X_test)
	test_sequences = [line.split('\t')[0] for line in open(args.test)]
	
	with open(prefix + '_predictions.txt', 'w') as outfile:
		for i in range(len(predictions)):
			outfile.write(
				test_sequences[i] + '\t' + 
				str(float(predictions[i])) + '\t' + 
				str(float(y_test[i])) + '\n')





