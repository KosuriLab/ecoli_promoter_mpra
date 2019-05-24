#!/usr/bin/python

from __future__ import division
import sys
from sklearn.neural_network import MLPRegressor
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


# def read_test_file(test_dataset_filename):
#     """
#     Load names and sequences for test data
#     """
#     sequences = []
#     names = []
#     with open(test_dataset_filename, 'r') as test_f:
#         for line in test_f:
#             data = (line.strip('\n')).split('\t')
#             sequence, expression = data[:]
#             names.append(name)
#             sequences.append(sequence)
#     return sequences, names


# def gen_dict(keys):
#     """
#     Generate dictionary given list of keys
#     """
#     kmer_dict = dict()
#     for key in keys:
#         kmer_dict[key] = 0
#     return kmer_dict


# def find_kmers(pos_kmers, string, k):
#     """
#     Count instances of kmers in sequence string
#     """
#     kmer_dict = gen_dict(pos_kmers)
#     str_length = len(string)
#     kmers = [string[i:i+k] for i in range(0, str_length-k+1)]
#     kmer_dict = dict(Counter(kmers))

#     return kmer_dict


# def find_rel_kmers(sequences, k, count_cutoff, use_cutoff):
#     """
#     Get kmer strings in training data meeting cutoff
#     """
#     kmer_dict = dict()
#     for seq in sequences:
#         seq_len = len(seq)
#         kmers = [seq[i:i+k] for i in range(0, seq_len-k+1)]
#         for kmer in kmers:
#             if kmer in kmer_dict:
#                 kmer_dict[kmer] += 1
#             else:
#                 kmer_dict[kmer] = 1
#     if use_cutoff:
#         cutoff = count_cutoff
#         print("cutoff", cutoff)
#         rel_kmer = []
#         for kmer in kmer_dict:
#             if kmer_dict[kmer] >= cutoff:
#                 rel_kmer.append(kmer)
#         return rel_kmer
#     else:
#         return sorted(kmer_dict.keys())


# def find_relevant_kmers(sequences, k, count_cutoff):
#     """
#     Get kmer strings in sequences above cutoff
#     """

#     kmer_counter = Counter()
#     for seq in sequences:
#         seq_len = len(seq)
#         kmers = [seq[i:i+k] for i in range(seq_len-k+1)]
#         kmer_counter.update(Counter(kmers))

#     rel_kmers = [x for x in kmer_counter if kmer_counter[x] >= count_cutoff]
#     return sorted(rel_kmers)


# def gen_features_kmers(sequences, k_min, k_max, test_sequences, count_cutoff):
    
#     """
#     Generate features from kmers
#     """

#     features = [[] for i in range(len(sequences))]
#     test_features = [[] for i in range(len(test_sequences))]
#     for k in range(k_min, k_max + 1):
#         relevant_kmers = find_relevant_kmers(sequences, k, count_cutoff)
#         print("k: ", k, ", len: ", len(relevant_kmers))
        
#         for idx, seq in enumerate(sequences):
#             kmer_dict = find_kmers(relevant_kmers, seq, k)
#             # for kmer in relevant_kmers:
#             # features[idx].append(kmer_dict[kmer])
#             sorted_kmer_counts = [kmer_dict.get(kmer, 0) for kmer in relevant_kmers]
#             features[idx] = sorted_kmer_counts

#         print "Created train features"
                
#         for idx, seq in enumerate(test_sequences):
#             kmer_dict = find_kmers(relevant_kmers, seq, k)
#             # for kmer in relevant_kmers:
#             #     test_features[idx].append(kmer_dict[kmer])
#             sorter_kmer_counts = [kmer_dict.get(kmer, 0) for kmer in relevant_kmers]
#             test_features[idx] = sorted_kmer_counts

#         print "Created test features"

#     return features, test_features


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



def predict(features, expression, test_features):
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
                # file_object = open(out_filename, "w")
                # line = "name\texpression\n"
                # file_object.write(line)
                # for idx, prediction in enumerate(predictions):
                #     line = '\t'.join([test_sequences[idx], str(prediction)]) + "\n"
                #     file_object.write(line)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('train', help='Filename of training set')
    parser.add_argument('test', help='Filename of test set')
    parser.add_argument('output_name', help='Filename of output')
    parser.add_argument('max_k', type=int, help='Max k-mer length')
    parser.add_argument('min_k', type=int, nargs='?', const=1, help='Minimum k-mer length')
    parser.add_argument('count_cutoff', type=int, nargs='?', const=0, help='Minimum k-mer count to count as feature')

    args = parser.parse_args()

    sequences, y_train = read_file(args.train)
    print 'loaded training data'
    test_sequences, y_test = read_file(args.test)
    print 'loaded test data'

X_train_all = pd.DataFrame()
X_test_all = pd.DataFrame()
# for k in range(args.min_k, args.max_k):
for k in range(min_k, max_k+1):
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
    scaler.fit(X_train)
    X_train_features = scaler.transform(X_train)
    X_test_features = scaler.transform(X_test)
    predictions = predict(X_train, y_train, X_test)

    print 'predicted expression'

    with open(args.output_name, 'w') as outfile:
        for i in range(len(predictions)):
            outfile.write(str(predictions[i]) + '\t' + str(y_test[i]) + '\n')

    print 'job complete'





