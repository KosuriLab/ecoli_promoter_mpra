from __future__ import absolute_import, division, print_function
import pandas as pd
import numpy as np, random
np.random.seed(1)
random.seed(1)
from dragonn.models import SequenceDNN
from dragonn.hyperparameter_search import HyperparameterSearcher, RandomSearch
try:
    from sklearn.model_selection import train_test_split  # sklearn >= 0.18
except ImportError:
    from sklearn.cross_validation import train_test_split  # sklearn < 0.18
from sklearn.preprocessing import LabelEncoder, OneHotEncoder
from sklearn.metrics import roc_curve, precision_recall_curve
import sys
import argparse
	

def one_hot_encode(sequences):
    sequence_length = len(sequences[0])
    integer_type = np.int8 if sys.version_info[
        0] == 2 else np.int32  # depends on Python version
    integer_array = LabelEncoder().fit(np.array(('ACGTN',)).view(integer_type)).transform(
        sequences.view(integer_type)).reshape(len(sequences), sequence_length)
    one_hot_encoding = OneHotEncoder(
        sparse=False, n_values=5, dtype=integer_type).fit_transform(integer_array)

    return one_hot_encoding.reshape(
        len(sequences), 1, sequence_length, 5).swapaxes(2, 3)[:, :, [0, 1, 2, 4], :]


def reverse_complement(encoded_seqs):
    return encoded_seqs[..., ::-1, ::-1]


def pad_sequence(seq, max_length):
	if len(seq) > max_length:
		diff = len(seq) - max_length
		# diff%2 returns 1 if odd
		trim_length = int(diff / 2)
		seq = seq[trim_length : -(trim_length + diff%2)]
	else:
		seq = seq.center(max_length, 'N')

	return seq

# def encode_trim_pad_fasta_sequences(fname, max_length):
#     """
#     One hot encodes sequences in fasta file. If sequences are too long, they will
#     be trimmed to the center. If too short, they will be padded with Ns
#     """
#     name, seq_chars = None, []
#     sequences = []
#     with open(fname) as fp:
#         for line in fp:
#             line = line.rstrip()
#             if line.startswith(">"):
#                 if name:
#                 	seq = ''.join(seq_chars).upper()
#                 	# this will center the string, and pad with Ns
#                 	if len(seq) > max_length:
#                 		diff = len(seq) - max_length
#                 		# diff%2 returns 1 if odd
#                 		trim_length = int(diff / 2)
#                 		seq = seq[trim_length : -(trim_length + diff%2)]
#                 	else:
#                 		seq = seq.center(max_length, 'N')
#                 	sequences.append(seq)
#                 name, seq_chars = line, []
#             else:
#                 seq_chars.append(line)
#     if name is not None:
#     	seq = ''.join(seq_chars).upper()
#     	# this will center the string, and pad with Ns
#     	if len(seq) > max_length:
#     		diff = len(seq) - max_length
#     		# diff%2 returns 1 if odd
#     		trim_length = int(diff / 2)
#     		seq = seq[trim_length : -(trim_length + diff%2)]
#     	else:
#     		seq = seq.center(max_length, 'N')
#         sequences.append(seq)

#     return one_hot_encode(np.array(sequences))


# def process_seqs(filename, seq_length, activity_type):

# 	X = encode_trim_pad_fasta_sequences(filename, seq_length)
# 	if activity_type == 'active':
# 		y = np.array([[True]]*len(X))
# 	elif activity_type == 'inactive':
# 		y = np.array([[False]]*len(X))
# 	else:
# 		raise ValueException('Please specify activity type: active or inactive')

# 	return [X, y]


def process_seqs(filename, seq_length):

	df = pd.read_csv(filename, header=None, names=['sequence', 'label'], sep='\t')
	padded_seqs = [pad_sequence(x, seq_length) for x in df['sequence'].tolist()]
	X = one_hot_encode(np.array(padded_seqs))
	y = np.array(df['label']).reshape(len(df), 1)
	# convert to bool
	y_bool = y == 1
	return [X, y_bool]


num_epochs = 100

if __name__ == '__main__':

	parser = argparse.ArgumentParser()
	parser.add_argument('train', help='''tab-separated, first column sequence, 
		second column label''')
	parser.add_argument('test', help='''tab-separated, first column sequence, 
		second column label''')

	parser.add_argument('seq_length', type=int, help='length of input sequences, or max length for trimming/padding sequences')
	parser.add_argument('num_layers', type=int, help='number of convolutional layers')
	parser.add_argument('min_filter', type=int, help='minimum number of filters')
	parser.add_argument('max_filter', type=int, help='maximum number of filters')
	parser.add_argument('validation_fraction', type=float)
	parser.add_argument('num_trials', type=int, 
		help='number of hyperparameter trials')
	parser.add_argument('prefix', help='output prefix for saved model files')
	args = parser.parse_args()

	train = args.train
	test = args.test
	seq_length = args.seq_length
	num_layers = args.num_layers
	min_filter = args.min_filter
	max_filter = args.max_filter
	validation_fraction = args.validation_fraction
	num_hyperparameter_trials = args.num_trials
	prefix = args.prefix

	# read in sequences and labels
	print("loading sequence data...")

	X_train, y_train = process_seqs(train, seq_length)

	X_train, X_valid, y_train, y_valid = train_test_split(X_train, y_train, 
		test_size=validation_fraction)

	X_test, y_test = process_seqs(test, seq_length)
	
	print('Starting hyperparameter search...')

	min_conv_width = 6
	max_conv_width = 30
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
	print('Test results: {}'.format(model.test(X_test, y_test)))
	
	# save model
	model.save(prefix)

	# predictions
	predictions = model.predict(X_test)

	with open(prefix + '_predictions.txt', 'w') as outfile:
		for i in range(len(predictions)):
			outfile.write(
				str(float(predictions[i])) + '\t' + 
				str(float(y_test[i])) + '\n')


	# fpr, tpr, thresholds = roc_curve(y_test, predictions)
	# with open(prefix + '_roc_info.txt', 'w') as outfile:
	# 	for i in range(len(fpr)):
	# 		outfile.write(str(fpr[i]) + ',' + str(tpr[i]) + ',' + str(thresholds[i]) + '\n')

	# precision, recall, thresholds = precision_recall_curve(y_test, predictions)
	# with open(prefix + '_pr_info.txt', 'w') as outfile:
	# 	for i in range(len(thresholds)):
	# 		outfile.write(str(precision[i]) + ',' + str(recall[i]) + ',' + str(thresholds[i]) + '\n')



