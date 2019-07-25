import argparse
from itertools import islice
import string
import numpy as np
import re
import random

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
		seqs[name] = seq
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


def generate_random_sequence(length, re_seqs):
	re_present = True

	while re_present:
		random_seq = ''.join([random.choice('ACGT') for i in range(length)])
		re_present = any([x in random_seq for x in re_sites])

	return random_seq


if __name__ == '__main__':
	parser = argparse.ArgumentParser("script to generate tiled regions of peak calls")
	parser.add_argument('sequences', help='filename of input sequences, FASTA format')
	parser.add_argument('neg_controls', help='fasta file of negative controls (more than 200bp from TSS either strand)')
	parser.add_argument('pos_controls', help='fasta file of positive controls')
	parser.add_argument('stride', type=int, help='distance between consecutive tiles')
	parser.add_argument('tile_len', type=int, help='tile length')
	parser.add_argument('output_name', help='name of output file')
	parser.add_argument('re_sites', help='fasta file, REs used for cloning, will remove sequences with RE sites')
	parser.add_argument('--n_random', type=int, help='Number of random sequences, optional')
	parser.add_argument('--rand_length', type=int, help='length of random sequences')

	args = parser.parse_args()
	sequences = fasta_reader(args.sequences)
	neg_controls = tab_reader(args.neg_controls)
	pos_controls = parse_controls(args.pos_controls)
	stride = args.stride
	tile_len = args.tile_len
	re_sites = fasta_reader(args.re_sites).values()
	if args.n_random:
		n_random = args.n_random
		rand_length = args.rand_length
	output_name = args.output_name

	# these primers include 20bp of primer and overlaps 2bp with 6bp restriction site
	# skpp-41-F sk20mer-106328
	fwd_primer = 'AAGACTCAACCAATGACCCTCGAG'
	# skpp-32-R sk20mer-25892
	rev_primer = 'GCTAGCCTATGTGTCCTACGAAGA'
	xhoI = 'CTCGAG'
	

	tiles = {}
	# T is easiest to synthesize, won't have secondary structure
	stuffer = 'T' * tile_len

	for name, seq in sequences.items():
		# if short, stuff
		if len(seq) < tile_len:
			# store as tuple so we can easily insert the restriction site
			# in between the stuffer and tile
			tile = add_stuffer(seq, stuffer, tile_len)
			tile_name = name + '_pos0-' + str(len(seq))
			tiles[tile_name] = tile
		else:
			for i in range(0, len(seq), stride):

				# if tile exceeds sequence length
				if i + tile_len > len(seq):
					start = len(seq) - tile_len
					end = len(seq)

				# # make sure full tile length contained in seq
				# if i + tile_len > (len(seq) - stride):
				# 	continue
				else:
					start = i
					end = i + tile_len

				tile = seq[start:end]
				tile_name = name + '_pos' + str(start) + '-' + str(end)
				tiles[tile_name] = tile

	# add controls to tiles
	tiles.update(neg_controls)
	# positive controls are shorter, need to stuff
	for name, seq in pos_controls.items():
		if len(seq) < tile_len:
			tiles[name] = add_stuffer(seq, stuffer, tile_len)
		else:
			tiles[name] = seq

	if args.n_random:
		# add random sequences
		random_seqs = {}
		for i in range(n_random):
			random_seq = generate_random_sequence(rand_length, re_sites)
			name = 'random' + str(i)
			tiles[name] = random_seq


	outfile = open(output_name, 'w')			
	# add restriction sites and primer
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

	outfile.close()





