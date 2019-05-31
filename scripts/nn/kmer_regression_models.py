from __future__ import division
import timeit
import sys
from sklearn.cross_decomposition import PLSRegression
from sklearn.preprocessing import StandardScaler
import argparse
from collections import Counter
import pandas as pd
import numpy as np
import statsmodels.api as sm


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('train', help='Filename of training set, k-mer counts')
    parser.add_argument('test', help='Filename of test set, k-mer counts')
    parser.add_argument('output_name', help='Output filename for predictions')
    parser.add_argument('method', help='regression method. linear, pls, svmlinear, svmradial')

    args = parser.parse_args()
    train_name = args.train
    test_name = args.test
    output_name = args.output_name
    method = args.method
    
    print "Reading training set..."
    X_train = pd.read_table(train_name)
    y_train = X_train['y']
    X_train.drop('y', axis=1, inplace=True)

    print "Reading test set..."
    X_test = pd.read_table(test_name)
    y_test = X_test['y']
    X_test.drop('y', axis=1, inplace=True)

    print "Fitting model..."
    if method == 'linear':
    	model = sm.OLS(y_train, X_train).fit()
    elif method == 'pls':
    	model = PLSRegression(n_components=2)
    	model.fit(X_train, y_train)

    print "Predicting..."
    predictions = model.predict(X_test)
    # check if list needs to be flattened
    if type(predictions[0]) != np.float64:
    	predictions = [item for sublist in predictions for item in sublist]
    score = np.corrcoef(y_test, predictions)[0, 1]
    print("Correlation:", score)

    print "Writing..."
    with open(output_name, 'w') as outfile:
    	for i in range(len(predictions)):
    		outfile.write(str(predictions[i]) + '\t' + str(y_test[i]) + '\n')

