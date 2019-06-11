#!/usr/bin/python

from __future__ import division
import timeit
import sys
from sklearn.neural_network import MLPRegressor, MLPClassifier
from sklearn.preprocessing import StandardScaler
import argparse
from collections import Counter
import pandas as pd


def read_file(dataset_filename):
    """
    Load sequences and expression values for training data
    """
    sequences = []
    expression = []
    with open(dataset_filename, 'r') as input_f:
        for line in input_f:
            data = (line.strip('\n')).split('\t')
            sequence, exp_value = data[:]
            sequences.append(sequence)
            exp = float(exp_value)
            expression.append(exp)
    return sequences, expression


def generate_filtered_kmer_counts(sequences, k, relevant_kmers=None, cutoff=None):
    
    # all_kmers = [''.join(x) for x in itertools.product(['A', 'C', 'G', 'T'], repeat=k)]
    
    start = timeit.default_timer()
    print "Counting k-mers of length", k, "in all sequences..."
    kmer_counts = [dict(Counter([seq[i:i+k] for i in range(len(seq)-k+1)])) for seq in sequences]
    elapsed = timeit.default_timer() - start
    print "Elapsed time:", elapsed
    # ensures all columns are present, even if some keys not present in all dictionaries
    print "Creating features..."
    kmer_df = pd.DataFrame.from_records(kmer_counts)
    if not relevant_kmers:
        kmer_sums = kmer_df.sum()
        relevant_kmers = list(kmer_sums[kmer_sums >= cutoff].index)
    print "Filtering..."
    # kmer_df_filtered = kmer_df.loc[:, kmer_df.columns.str.contains('|'.join(relevant_kmers))]
    kmer_df_filtered = kmer_df[relevant_kmers]
    kmer_df_filtered.fillna(0, inplace=True)
    print "Number of filtered k-mers:", len(relevant_kmers)
    return kmer_df_filtered, relevant_kmers



def predict_regression(features, expression, test_features):
    """
    Predict expression of test sequences
    """
    for a in [5e-3]:
        for input_layer in [800]:
            for hidden_layer in [30]:
                print("alpha: ", a, 
                    " input_layer: ", input_layer, 
                    " hidden_layer: ", hidden_layer)
                clf = MLPRegressor(solver="lbfgs", 
                    alpha=a, 
                    hidden_layer_sizes=(input_layer, hidden_layer), 
                    random_state=1, 
                    max_iter=10000, 
                    verbose=False, 
                    early_stopping=True, 
                    learning_rate="adaptive", 
                    tol=1e-8)

                clf.fit(features, expression)
                predictions = clf.predict(test_features)

    return predictions


def predict_classification(features, expression, test_features):
    """
    Predict expression of test sequences
    """
    for a in [5e-3]:
        for input_layer in [800]:
            for hidden_layer in [30]:
                print("alpha: ", a, 
                    " input_layer: ", input_layer, 
                    " hidden_layer: ", hidden_layer)
                clf = MLPClassifier(solver="lbfgs", 
                    alpha=a, 
                    hidden_layer_sizes=(input_layer, hidden_layer), 
                    random_state=1, 
                    max_iter=10000, 
                    verbose=False, 
                    early_stopping=True, 
                    learning_rate="adaptive", 
                    tol=1e-8)

                clf.fit(features, expression)
                predictions = clf.predict(test_features)

    return predictions


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('train', help='Filename of training set')
    parser.add_argument('test', help='Filename of test set')
    parser.add_argument('output_name', help='Filename of output')
    parser.add_argument('min_k', type=int, help='Minimum k-mer length')
    parser.add_argument('max_k', type=int, help='Maximum k-mer length')
    parser.add_argument('count_cutoff', type=int, nargs='?', const=0, help='Minimum k-mer count to count as feature')
    parser.add_argument('--regression', action='store_true')
    parser.add_argument('--classification', action='store_true')

    args = parser.parse_args()
    
    count_cutoff = args.count_cutoff
    sequences, y_train = read_file(args.train)
    print 'loaded training data'
    test_sequences, y_test = read_file(args.test)
    print 'loaded test data'

    X_train_all = pd.DataFrame()
    X_test_all = pd.DataFrame()
    for k in range(args.min_k, args.max_k+1):
    # for k in range(min_k, max_k+1):
        print "Train:"
        X_train, relevant_kmers = generate_filtered_kmer_counts(sequences, k, cutoff=count_cutoff)
        X_train_all = pd.concat([X_train_all, X_train], axis=1)
        print "Test:"
        X_test, relevant_kmers = generate_filtered_kmer_counts(test_sequences, k,
            relevant_kmers=relevant_kmers)
        X_test_all = pd.concat([X_test_all, X_test], axis=1)

    print 'generated features'
    # remove mean and scale to unit variance
    scaler = StandardScaler()
    scaler.fit(X_train_all)
    X_train_features = scaler.transform(X_train_all)
    X_test_features = scaler.transform(X_test_all)
    if args.regression:
        predictions = predict_regression(X_train_features, y_train, X_test_features)
    if args.classification:
        predictions = predict_classification(X_train_features, y_train, X_test_features)

    print 'predicted expression'

    with open(args.output_name, 'w') as outfile:
        for i in range(len(predictions)):
            outfile.write(str(predictions[i]) + '\t' + str(y_test[i]) + '\n')

    print 'job complete'





