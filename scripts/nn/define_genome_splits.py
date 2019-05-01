import argparse
import random
import pandas as pd
from math import floor
from copy import deepcopy


def create_splits(train_fraction):

	test_fraction = 1 - train_fraction
	split_fraction = test_fraction / 2
	split_size = int(floor(genome_size * split_fraction))

	# genome coordinates are 1-based
	splits = [(i, i + split_size - 1) for i in range(1, genome_size + 1, split_size)]

	# remove last split which contains the remainder of the split size
	splits.pop()
	# fix last remaining split
	last_split = splits.pop()
	splits.append((last_split[0], genome_size))

	return splits


def find_equidistant_pair(splits):
	'''
	Input: list of tuples corresponding to the start and end of the genome split
	Output: two splits that are equidistant from the origin. The genome is circular
	and the origin is at coordinate 1 and coordinate <genome size>, the second 
	element in the last tuple
	'''


	start = splits[0][0]
	end = splits[-1][1]

	# to accurately compare distances, determine how much longer the last split
	# is since it contains the modulus of genome % split_size
	split_size = splits[0][1] - splits[0][0] + 1
	remainder = end % split_size

	# randomly choose split
	split1 = random.choice(splits)
	# determine distance from origin
	distance1 = split1[0] - start 
	# iterate through splits until we find out that is equidistant
	for i in range(len(splits)):
		split2 = splits[i]
		distance2 = abs(split2[1] - end + remainder)
		if distance1 == distance2:
			return [split1, split2]

	return None


def in_range(split,  x, y):
	# Check if (x, y) is contained in split, a tuple (start, end)
	if x >= split[0] and y <= split[1]:
		return True
	else:
		return False


if __name__ == '__main__':

	parser = argparse.ArgumentParser('''Define genome positions used for train/test.
		Test and train locations will be equidistant from origin of replication.''')
	parser.add_argument('train_fraction', type=float, help='''Fraction of genomic
		positions used for training''')
	parser.add_argument('genome_size', type=int, help='Size of genome (bp)')
	parser.add_argument('infile_train',  help='''filename of formatted scramble expression data.
		Must include columns for variant, expn_med_fitted_scaled, var_left, var_right''')
	parser.add_argument('infile_test', help='''filename of formatted TSS expression data.
		Must include columns for variant, expn_med_fitted_scaled, start, end''')
	parser.add_argument('outfile_train', help='''Output name for position split
		train set. Two-column, tab-separated''')
	parser.add_argument('outfile_test', help='''Output name for position split
		test set. Two-column, tab-separated''')
	args = parser.parse_args()

	train_fraction = args.train_fraction
	genome_size = args.genome_size

	splits = create_splits(train_fraction)	

	# assign test splits so they are equidistant from the origin (coordinate 1, circular)
	# we created the splits so two splits are equal to the test fraction
	# we must choose two splits equidistant and assign to test, the rest to train
	test_splits = find_equidistant_pair(splits)
	# create train/test look-up for each split
	split_lookup = {x : 'test' if x in test_splits else 'train' for x in splits}

	# read in datasets
	train = pd.read_table(args.infile_train, sep='\t', header=0)
	test = pd.read_table(args.infile_test, sep='\t', header=0)

	# sort by genomic position, make sure all positive
	train.sort_values(by='var_left', inplace=True)
	test.sort_values(by='start', inplace=True)
	test = test[test.start > 0]

	# keep relevant categories
	train = train[train.category == 'scramble']
	test = test[test.category == 'tss']

	# For train, only sequences that fall into train splits will be written
	with open(args.outfile_train, 'w') as train_outfile:
	# with open("train_test.txt", 'w') as train_outfile:

		for i in range(len(train)):
			x = train.var_left.iloc[i]
			y = train.var_right.iloc[i]
			
			for j in range(len(splits)):
				if in_range(splits[j], x, y):
					if split_lookup[splits[j]] == 'train':
						train_outfile.write(train.variant.iloc[i] + '\t')
						train_outfile.write(str(train.expn_med_fitted_scaled.iloc[i]) + '\n')
					# else, don't write to output


	# For test, only sequences that fall into test splits will be written
	with open(args.outfile_test, 'w') as test_outfile:
	# with open("test_test.txt", 'w') as test_outfile:
		for i in range(len(test)):
			x = test.start.iloc[i]
			y = test.end.iloc[i]
			for j in range(len(splits)):
				if in_range(splits[j], x, y):
					if split_lookup[splits[j]] == 'test':
						test_outfile.write(test.variant.iloc[i] + '\t')
						test_outfile.write(str(test.expn_med_fitted_scaled.iloc[i]) + '\n')
					# else, don't write to output



