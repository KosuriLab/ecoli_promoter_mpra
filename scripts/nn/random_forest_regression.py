# from __future__ import absolute_import, division, print_function
import pandas as pd
import numpy as np, random
np.random.seed(1)
random.seed(1)
from keras_regression import SequenceDNN, RandomForestRegression, DecisionTree
from hyperparameter_search_regression import HyperparameterSearcher, RandomSearch
from sklearn.preprocessing import LabelEncoder, OneHotEncoder
try:
    from sklearn.model_selection import train_test_split  # sklearn >= 0.18
except ImportError:
    from sklearn.cross_validation import train_test_split  # sklearn < 0.18
import sys
import argparse
import time
import itertools
from collections import Counter

num_epochs = 100



def process_seqs(filename):

	data = pd.read_table(filename, header=None)
	sequences = data.loc[:, 0]
	y = np.array(data.loc[:, 1])
	X = np.array(data.loc[:, 2:])
	
	return X, y





if __name__ == '__main__':

	parser = argparse.ArgumentParser()
	parser.add_argument('train', help='''Training set, one column sequence, 
		one column expression. Tab-separated''')
	parser.add_argument('test', help='''Test sest, one column sequence, 
		one column expression. Tab-separated''')
	parser.add_argument('output_name')
	args = parser.parse_args()

	# load in pre-defined splits
	print("loading training and test set...")
	X_train, y_train = process_seqs(args.train)
	X_test, y_test = process_seqs(args.test)

	print("Running random forest regression...")
	model = RandomForestRegression()
	model.train(X_train, y_train)
	predictions = model.predict(X_test)

	with open(args.output_name, 'w') as outfile:
		for i in range(len(predictions)):
			outfile.write(str(float(predictions[i])) + '\t' +
				      str(float(y_test[i])) + '\n')

	score = model.score(X_test, y_test)
	print("Score:", score)

