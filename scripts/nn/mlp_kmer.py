#!/usr/bin/python

from __future__ import division
import sys
from sklearn.neural_network import MLPRegressor
from sklearn.preprocessing import StandardScaler
import argparse

def read_file(dataset_filename):
    """
    Load sequences and expression values for training data
    """
    sequences = []
    expression = []
    with open(dataset_filename, 'r') as input_f:
        for line in input_f:
            data = (line.strip('\n')).split('\t')
            name, sequence, exp_value = data[:]
            sequences.append(sequence)
            exp = float(exp_value)
            expression.append(exp)
    return sequences, expression


def read_test_file(test_dataset_filename):
    """
    Load names and sequences for test data
    """
    sequences = []
    names = []
    with open(test_dataset_filename, 'r') as test_f:
        for line in test_f:
            data = (line.strip('\n')).split('\t')
            name, sequence, expression = data[:]
            names.append(name)
            sequences.append(sequence)
    return sequences, names


def gen_dict(keys):
    """
    Generate dictionary given list of keys
    """
    kmer_dict = dict()
    for key in keys:
        kmer_dict[key] = 0
    return kmer_dict


def find_kmers(pos_kmers, string, k):
    """
    Count instances of kmers in sequence string
    """
    kmer_dict = gen_dict(pos_kmers)
    str_length = len(string)
    kmers = [string[i:i+k] for i in range(0, str_length-k+1)]
    for kmer in kmers:
        if kmer in kmer_dict:
            kmer_dict[kmer] += 1
    return kmer_dict


def find_rel_kmers(sequences, k, count_cutoff, use_cutoff):
    """
    Get kmer strings in training data meeting cutoff
    """
    kmer_dict = dict()
    for seq in sequences:
        seq_len = len(seq)
        kmers = [seq[i:i+k] for i in range(0, seq_len-k+1)]
        for kmer in kmers:
            if kmer in kmer_dict:
                kmer_dict[kmer] += 1
            else:
                kmer_dict[kmer] = 1
    if use_cutoff:
        cutoff = count_cutoff
        print("cutoff", cutoff)
        rel_kmer = []
        for kmer in kmer_dict:
            if kmer_dict[kmer] >= cutoff:
                rel_kmer.append(kmer)
        return rel_kmer
    else:
        return sorted(kmer_dict.keys())


def gen_features_kmers(sequences, k_min, k_max, test_sequences, count_cutoff, use_cutoff):
    """
    Generate features from kmers
    """
    features = [[] for i in range(len(sequences))]
    test_features = [[] for i in range(len(test_sequences))]
    for k in range(k_min, k_max + 1):
        pos_kmers = find_rel_kmers(sequences, k, count_cutoff, use_cutoff)
        print("k: ", k, ", len: ", len(pos_kmers))
        for idx, seq in enumerate(sequences):
            kmer_dict = find_kmers(pos_kmers, seq, k)
            for kmer in pos_kmers:
                features[idx].append(kmer_dict[kmer])
        for idx, seq in enumerate(test_sequences):
            kmer_dict = find_kmers(pos_kmers, seq, k)
            for kmer in pos_kmers:
                test_features[idx].append(kmer_dict[kmer])
    return features, test_features


def predict(features, expression, test_features, test_names, out_filename, test_sequences):
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
                file_object = open(out_filename, "w")
                line = "name\texpression\n"
                file_object.write(line)
                for idx, prediction in enumerate(predictions):
                    line = '\t'.join([test_names[idx], test_sequences[idx], str(prediction)]) + "\n"
                    file_object.write(line)


def main(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument('train', help='Filename of training set')
    parser.add_argument('test', help='Filename of test set')
    parser.add_argument('output_name', help='Filename of output')
    parser.add_argument('max_k', type=int, help='Max k-mer length')
    parser.add_argument('min_k', type=int, const=1, help='Minimum k-mer length')
    parser.add_argument('count_cutoff', type=int, const=100, help='Minimum k-mer count to count as feature')

    args = parser.parse_args()

    sequences, expression = read_file(args.train)
    print 'loaded training data'
    test_sequences, names = read_test_file(args.test)
    print 'loaded test data'

    features, test_features = gen_features_kmers(
        sequences=sequences, 
        k_min=args.min_k, 
        k_max=args.max_k,
        test_sequences=test_sequences, 
        count_cutoff=args.count_cutoff, 
        use_cutoff=True)
    print 'generated features'

    scaler = StandardScaler()
    scaler.fit(features)
    features = scaler.transform(features)
    test_features = scaler.transform(test_features)
    predict(features, expression, test_features, names, args.output_name, test_sequences)
    print 'predicted expression'

    print 'job complete'



if __name__ == "__main__":
    main(sys.argv[1:])



