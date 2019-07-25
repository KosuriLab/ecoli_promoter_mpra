import argparse
from itertools import islice
import string
import numpy as np
import random
import re
from difflib import SequenceMatcher

random.seed(123)

def fasta_reader(filename):
	"""
	Input: str, name of file
	Output: dictionary, key = header, value = sequence
	"""

	seqs = {}
	with open(filename) as infile:
		entry = list(islice(infile, 2))
		while entry:
			# grab two lines from file at a time, strip \n
			header, seq = map(string.strip, entry)
			# strip '>' from header, first character
			seqs[header[1:]] = seq.upper()
			entry = list(islice(infile, 2))

		return seqs


def tab_reader(filename):
	infile = open(filename)
	seqs = {}

	for line in infile:
		name, seq = line.strip().split('\t')
		# remove leading '>' from fasta-style header if present
		name = name.replace('>', '')
		seqs[name] = seq

	return seqs


def parse_controls(filename):
	'''
	Extract promoters from CSV file. There are weird quotation marks and spaces
	so have to do some regex stuff
	'''

	infile = open(filename, 'r')

	syn_promoters = {}

	# read through header
	infile.readline()

	for line in infile.readlines():
		fields = line.strip().split(',')
		name = fields[0]
		seq = fields[9]

		# remove white space from seq
		seq = ''.join(seq.split())
		# remove quotations
		match = re.search('[ACGT]{1,}', seq)
		clean_seq = match.group(0)

		# remove quotations from name, always first and last characters
		clean_name = name.replace('\"', '')

		syn_promoters['pos_control_'+clean_name] = clean_seq

	return syn_promoters


def best_A_content(oligo):
	'''
	Choose the strand with the lowest A content because A's are harder to
	synthesize. Return true if sequence needs to be reverse complemented
	'''
	rc_oligo = reverse_complement(oligo)

	oligo_As = sum( [1 for nt in oligo if nt == 'A'] )
	rc_As = sum( [1 for nt in rc_oligo if nt == 'A'] )

	if oligo_As < rc_As:
		return False
	else:
		return True


def reverse_complement(sequence):
    """Return the reverse complement of the supplied sequence string """ 
    basecomplement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N':'N', '.':'.', '*':'*'} 
    #Reverse, convert to uppercase
    sequence = sequence[::-1].upper()
    #Complement
    letters = list(sequence) 
    letters = [basecomplement[base] for base in letters] 
    return ''.join(letters) 


def add_stuffer(sequence, stuffer, tile_len):
	# store as tuple so we can easily insert the restriction site in between the
	# stuffer and tile
	stuffer = (stuffer[:(tile_len - len(sequence))], sequence)
	return stuffer


def similar(a, b):
    # must be same length sequence
    num_same = sum([1 for i in range(len(a)) if a[i] == b[i]])
    return num_same


def most_scramble(sequence, n_iter):
	'''
	Scramble sequence n_iter times and choose the most scrambled sequence
	'''
	scrambled_seqs = [''.join(random.sample(list(sequence), len(sequence))) for i in range(n_iter)]
	similarity = [similar(sequence, x) for x in scrambled_seqs]
	max_scrambled = scrambled_seqs[np.argmin(similarity)]
	return max_scrambled


if __name__ == '__main__':
	parser = argparse.ArgumentParser("script to generate scrambled regions tiling input sequences")
	parser.add_argument('sequences', help='filename of input sequences, FASTA format')
	parser.add_argument('neg_controls', help='fasta file of negative controls (more than 200bp from TSS either strand)')
	parser.add_argument('pos_controls', help='fasta file of positive controls')
	parser.add_argument('scramble_len', type=int, help='length of scrambled segments')
	parser.add_argument('stride_len', type=int, help='step length between scrambled segments')
	
	parser.add_argument('output_name', help='name of output file')

	args = parser.parse_args()
	sequences = fasta_reader(args.sequences)
	neg_controls = tab_reader(args.neg_controls)
	pos_controls = parse_controls(args.pos_controls)
	scramble_len = args.scramble_len
	stride_len = args.stride_len
	output_name = args.output_name

	# just grab length of input sequence
	tile_len = len(sequences.values()[0])
	# T is easiest to synthesize, won't have secondary structure
	stuffer = 'T' * tile_len

	# skpp-100-F sk20mer-20955
	fwd_primer = 'ACCTGTAATTCCAAGCGTCTCGAG'
	# skpp-155-R sk20mer-121927
	rev_primer = 'GCTAGCGGTGTTTAGTTAGCATCC'
	xhoI = 'CTCGAG'

	tiles = {}

	for name, seq in sequences.items():
		# add unscrambled
		tile_name = name + '_unscrambled'
		tiles[tile_name] = seq
		for i in range(0, len(seq), stride_len):
			if i + scramble_len > len(seq):
				continue
			# scrambled = list(seq[i:i+scramble_len])
			# random.shuffle(scrambled)
			to_scramble = seq[i:i+scramble_len] 
			scrambled = most_scramble(to_scramble, 100)
			tile = seq[:i] + ''.join(scrambled).lower() + seq[i+scramble_len:]
			# print ''.join(scrambled).lower()
			tile_name = name + '_scrambled' + str(i) + '-' + str(i+scramble_len)
			tiles[tile_name] = tile

	# add controls to tiles
	tiles.update(neg_controls)
	# positive controls are shorter, need to stuff
	for name, seq in pos_controls.items():
		if len(seq) < tile_len:
			tiles[name] = add_stuffer(seq, stuffer, tile_len)
		else:
			tiles[name] = seq


	with open(output_name, 'w') as outfile:
		for tile_name in sorted(tiles.keys()):
			tile = tiles[tile_name]

			if len(tile) == 2:
				stuffer, payload = tile
				# stuffed sequence, only use first 20bp of fwd primer, subtract 2bp from stuffer,
				# then add the full 6bp of RE (no overlap between primer and RE anymore)
				full_tile = fwd_primer[:20] + stuffer[:-2] + xhoI + payload + rev_primer
			else:
				full_tile = fwd_primer + tile + rev_primer

				# 149bp peaks end up being 1bp too long, remove bp from rev primer
			if len(full_tile) == 199:
				full_tile = fwd_primer[:20] + stuffer[:-2] + xhoI + payload + rev_primer[:23]
			
			reverse = best_A_content(full_tile)
			if reverse:
				tile_name += '_flipped'
				full_tile = reverse_complement(full_tile)

			outfile.write(tile_name + '\t' + full_tile + '\n')


