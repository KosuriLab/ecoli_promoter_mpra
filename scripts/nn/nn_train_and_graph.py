from __future__ import absolute_import, division, print_function
import numpy as np, random
np.random.seed(1)
random.seed(1)
from keras_regression import SequenceDNN, SequenceDNN_dropout
from hyperparameter_search_regression import HyperparameterSearcher, RandomSearch
from sklearn.preprocessing import LabelEncoder, OneHotEncoder
try:
    from sklearn.model_selection import train_test_split  # sklearn >= 0.18
except ImportError:
    from sklearn.cross_validation import train_test_split  # sklearn < 0.18
import sys
import argparse
import matplotlib.pyplot as plt


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


def encode_trim_pad_fasta_sequences(fname, max_length):
    """
    One hot encodes sequences in fasta file. If sequences are too long, they will
    be trimmed to the center. If too short, they will be padded with Ns
    """
    name, seq_chars = None, []
    sequences = []
    with open(fname) as fp:
        for line in fp:
            line = line.rstrip()
            if line.startswith(">"):
                if name:
                	seq = ''.join(seq_chars).upper()
                	# this will center the string, and pad with Ns
                	if len(seq) > max_length:
                		diff = len(seq) - max_length
                		# diff%2 returns 1 if odd
                		trim_length = int(diff / 2)
                		seq = seq[trim_length : -(trim_length + diff%2)]
                	else:
                		seq = seq.center(max_length, 'N')
                	sequences.append(seq)
                name, seq_chars = line, []
            else:
                seq_chars.append(line)
    if name is not None:
    	seq = ''.join(seq_chars).upper()
    	# this will center the string, and pad with Ns
    	if len(seq) > max_length:
    		diff = len(seq) - max_length
    		# diff%2 returns 1 if odd
    		trim_length = int(diff / 2)
    		seq = seq[trim_length : -(trim_length + diff%2)]
    	else:
    		seq = seq.center(max_length, 'N')
        sequences.append(seq)

    return one_hot_encode(np.array(sequences))


def process_seqs(filename, seq_length):

	seqs = [line.split('\t')[0] for line in open(filename)]
	padded_seqs = [pad_sequence(x, seq_length) for x in seqs]
	X = one_hot_encode(np.array(padded_seqs))
	y = np.array([float(line.strip().split('\t')[1]) for line in open(filename)])
	return X, y

	
if __name__ == '__main__':
	parser = argparse.ArgumentParser('Train keras sequential CNN with regression and plot correlation')
	parser.add_argument('sequences', help='tab-separated, two columns. First is sequence, second is continuous value')
	parser.add_argument('seq_length', type=int, help='length of input sequences')
	parser.add_argument('num_layers', type=int, help='number of convolutional layers')
	parser.add_argument('dropout', type=float, help='dropout rate')
	parser.add_argument('pool_width', type=int, help='width of max pooling layer')
	parser.add_argument('conv_width', type=str, 
		help='width of convolutional filter, comma separated string, one value for each layer')
	parser.add_argument('num_filters', type=str, 
		help='number of filters in each layer, comma separated string, one value for each layer')
	parser.add_argument('test_fraction', type=float, default=0.2, help='fraction used for testing')
	parser.add_argument('validation_fraction', type=float, default=0.2, help='fraction used for valdation')
	parser.add_argument('output_prefix', help='prefix of output graph')
	parser.add_argument('--train', help='pre-defined training set. Test fraction will be ignored but\
		validation fraction still applies to training set.')
	parser.add_argument('--test', help='pre-defined test set')
	parser.add_argument('--uncertain', action='store_true', default=False, help='predict with uncertainty')
	args = parser.parse_args()
 	
 	sequences = args.sequences
 	seq_length = args.seq_length
	num_layers = args.num_layers
	dropout = args.dropout
	pool_width = args.pool_width
	conv_width = map(int, args.conv_width.split(','))
	num_filters = map(int, args.num_filters.split(','))
	output_name = args.output_prefix
	test_fraction = args.test_fraction
	validation_fraction = args.validation_fraction
	uncertain = args.uncertain


	if args.train:
		# check test exists
		if not args.test:
			raise Exception("Please provide test set.")
		else:
			# load in pre-defined splits
			print("loading training set...")
			X_train, y_train = process_seqs(args.train, seq_length)
			print("loading test set...")
			X_test, y_test = process_seqs(args.test, seq_length)

	else:
		# read in sequences and labels
		print("loading sequence data...")
		X, y = process_seqs(sequences, seq_length)

		# need test index so we can grab these sequences later and output them for 
		# downstream analysis
		seqs = [line.split('\t')[0] for line in open(sequences)]
		random_test_index = random.sample(
			list(range(len(seqs))), 
			int(round(test_fraction * len(seqs)))
			)
		random_test_index = sorted(random_test_index)
		train_index = [i for i in range(len(seqs)) if i not in random_test_index]

		X_test = np.take(X, random_test_index, axis=0)
		y_test = np.take(y, random_test_index, axis=0)
		X_train = np.take(X, train_index, axis=0)
		y_train = np.take(y, train_index, axis=0)

	# split training into validation set
	X_train, X_valid, y_train, y_valid = train_test_split(X_train, y_train, 
		test_size=validation_fraction)

	# print("loading sequence data...")
	# seqs = [line.split('\t')[0] for line in open(sequences)]
	# padded_seqs = [pad_sequence(x, seq_length) for x in seqs]
	# X = one_hot_encode(np.array(padded_seqs))
	# # X = encode_trim_pad_fasta_sequences(sequences, seq_length)
	# y = np.array([float(line.strip().split('\t')[1]) for line in open(sequences)])

	# # need test index so we can grab these sequences later and output them for 
	# # downstream analysis
	# random_test_index = random.sample(
	# 	list(range(len(seqs))), 
	# 	int(round(test_fraction * len(seqs)))
	# 	)
	# random_test_index = sorted(random_test_index)
	# train_index = [i for i in range(len(seqs)) if i not in random_test_index]

	# X_test = np.take(X, random_test_index, axis=0)
	# y_test = np.take(y, random_test_index, axis=0)
	# X_train = np.take(X, train_index, axis=0)
	# y_train = np.take(y, train_index, axis=0)


	if uncertain:
		model = SequenceDNN_dropout(
			seq_length = seq_length,
			num_filters=num_filters,
			conv_width=conv_width,
			pool_width=pool_width,
			dropout=dropout)
	else:
		model = SequenceDNN(
			seq_length = seq_length,
			num_filters=num_filters,
			conv_width=conv_width,
			pool_width=pool_width,
			dropout=dropout)

	model.train(X_train, y_train, validation_data=(X_valid, y_valid))

	if uncertain:
		predictions, uncertainty = model.predict_with_uncertainty(
			X=X_test, 
			num_classes=1,
			n_iter=100,
			length_scale=1,
			dropout=dropout,
			L1=0)
	else:
		predictions = model.predict(X_test)


	corr = np.corrcoef(np.squeeze(predictions), y_test)[0, 1]
	print('Test results: {}'.format(corr))
	model.save(output_name + '_trained_model')

	if args.test:
		test_sequences = [line.split('\t')[0] for line in open(args.test)]
	else:
		test_sequences = [seqs[i] for i in random_test_index]
	
	if uncertain:
		with open(output_name + '_predictions.txt', 'w') as outfile:
			for i in range(len(predictions)):
				outfile.write(
					test_sequences[i] + '\t' + 
					str(float(predictions[i])) + '\t' + 
					str(float(y_test[i])) + '\t' +
					str(float(uncertainty[i])) + '\n')
	else:
		with open(output_name + '_predictions.txt', 'w') as outfile:
			for i in range(len(predictions)):
				outfile.write(
					test_sequences[i] + '\t' + 
					str(float(predictions[i])) + '\t' + 
					str(float(y_test[i])) + '\n')


	# predictions = np.squeeze(model.predict(X_test))
	# corr_text = 'r = ' + str(round(corr, 3))
	# plt.figure()
	# plt.scatter(y_test, predictions, alpha=0.50)
	# plt.text(0.5, 75, corr_text)
	# plt.yscale('symlog')
	# plt.xscale('symlog')
	# plt.title('Neural net predictions for held-out test set (n = ' + str(len(y_test)) + ')')
	# plt.xlabel('observed')
	# plt.ylabel('predicted')
	# plt.savefig(output_name + '.png')





