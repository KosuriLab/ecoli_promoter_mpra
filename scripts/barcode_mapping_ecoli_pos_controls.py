'''
This script separately maps the synthetic positive controls. They are
variable length and flanked by stuffer sequence, so it's not easy to extract
a fixed position with normal mapping scripts.
'''

import sys
import os
import itertools
from collections import Counter, defaultdict
import re
import Levenshtein
import numpy
import random
import argparse
import subprocess
import shlex


def fastq_reader(filename):
	"""
	Extract sequences from FASTQ file. Each sequence takes up four lines,
	with the sequence being the second line.
	"""

	with open(filename) as infile:
		# grab lines starting at second line (1, 0-indexing), go until end of
		# file (stop = None), skip 4 lines at a time 
		line_slice = itertools.islice(infile, 1, None, 4)
		
		# convert iterator to list, strip new lines
		sequences = []
		for line in line_slice:
			sequences.append(line.strip())

	return sequences


def fasta_reader(filename):

	seqs = {}

	with open(filename) as infile:
		seq = ''
		for line in infile:
			if line.startswith('>'):
				if len(seq) != 0:
					# print name, seq
					seqs[name] = seq
				
				name = line.strip()[1:] # remove leading '>'
				seq = ''
			else:
				seq += line.strip()

		# catch last sequence
		if len(seq) != 0:
			# print name, seq
			seqs[name] = seq

	return seqs


def sam_parser(filename):
	infile  = open(filename)
	aligned_reads = []
	alignments = {}
	for line in infile:
		fields = line.strip().split('\t')
		if fields[1] != '4': # read is mapped
			read = fields[9]
			name = fields[2]
			# # don't know why it reports 'null' even though
			# # > 50bp map to a reference sequence, just parse separately
			# if name == 'null':
			# 	# check cigar string to grab portion that matches, see if it
			# 	# matches
			# 	cigar = fields[5]
			# 	if cigar[2] == '=':
			# 		match_len = int(cigar[:2])
			# 		variant = read[:match_len]
			# 		if variant in controls:
			# 			# assign alignment proper name
			# 			alignments[read] = controls[variant]
			# 			aligned_reads.append(read)
			# else:
				
			alignments[read] = name
			aligned_reads.append(read)

	return aligned_reads, alignments



def library_reader(filename, primer_len, rev_complement=True, format = 'csv'):
	"""
	Read in .csv or tab file of library sequences. First column is sequence name,
	second column is sequence. Trim primer sequences.
	"""

	lib = {}

	with open(filename) as infile:
		# # read through first line
		# infile.readline()
		for line in infile:
			if format == 'csv':
				name, seq = line.strip().split(',')[:2]
			elif format == 'tab':
				name, seq = line.strip().split('\t')[:2]
			
			seq = seq.upper()
			if primer_len != 0:
				seq = seq[primer_len:-primer_len]

			if rev_complement:
				rc_seq = reverse_complement(seq)
				lib[rc_seq] = name + '_rc'

			lib[seq] = name

	return lib


def reverse_complement(seq):
	"""
	Return the reverse complement of a nucleotide string
	"""
	complement = {'A': 'T', 'T':'A', 'C':'G', 'G':'C'}
	
	rc = ''.join([complement[nt] for nt in seq[::-1]])
	return rc


def mapping(barcodes, aligned_reads, reads, alignments, controls, bc_loc, bc_len):

	variant_map = defaultdict(list)
	barcode_map = defaultdict(list)

	# create reverse dictionary of controls where names are keys
	controls_flipped = {controls[x] : x for x in controls}

	# for each barcode that passes filters, look at all reads and see what
	# it maps to
	for read in reads:
		if bc_loc == 'start':
			barcode = read[:bc_len]
		elif bc_loc == 'end':
			barcode = read[-bc_len:]

		# if barcode in filtered set
		if barcodes.get(barcode, 0) > 0:
			barcode_map[barcode].append(read)

	# in the aligned reads, keep track of how many barcodes go to each variant
	for read in aligned_reads:
		# look up aligned control
		name = alignments[read]
		# find proper variant length
		var_len = len(controls_flipped[name])
		if bc_loc == 'start':
			barcode = read[:bc_len]
			variant = read[bc_len : bc_len+var_len]
		elif bc_loc == 'end':
			variant = read[:var_len]
			barcode = read[-bc_len:]

		variant_map[variant].append(barcode)

	return [variant_map, barcode_map]


def filter_barcodes(barcode_map, controls, alignments, cutoff, bc_loc, bc_len, name):
	'''
	For each barcode, calculate the Levenshtein distance between its reads
	and if it is below cutoff (aka barcode maps to similar reads, no cross-talk)
	then keep this barcode
	'''

	final_barcodes = []
	covered_sequences = set()
	# controls = set(controls.keys()) # for faster lookup
	# create reverse dictionary of controls where names are keys
	controls_flipped = {controls[x] : x for x in controls}

	all_dist = []

	outfile = open(name, 'w')
	headers = ['barcode', 'num_unique', 'num_reads', 'num_reads_most_common', 
	'most_common', 'name']
	outfile.write('\t'.join(headers)+'\n')

	for barcode in barcode_map:
		reads = barcode_map[barcode]
		# grab most common read as reference
		most_common = Counter(reads).most_common(1)[0][0]

		# find positive control, shorten trimmed reads to appropriate length
		name = alignments.get(most_common, None)
		if name is not None:
			var_len = len(controls_flipped[name])
			trimmed = [read[:var_len] for read in reads]
			most_common_trimmed = most_common[:var_len]

			# calculate distance between each read and reference
			distances = [Levenshtein.distance(most_common_trimmed, read) for read in set(trimmed)]
			all_dist.append(max(distances))
			# max distance for set of reads belonging to a barcode must be below cutoff
			if max(distances) < cutoff:
				num_unique = len(set(trimmed))
				num_reads = len(trimmed)
				num_reads_most_common = Counter(trimmed).most_common(1)[0][1]
				covered_sequences.add(name)
				final_barcodes.append(barcode)
				info = [barcode, num_unique, num_reads, num_reads_most_common, most_common_trimmed, name]
				info = map(str, info)
				outfile.write('\t'.join(info)+'\n')

	outfile.close()

	coverage = len(covered_sequences)/ float((len(controls)/2.0))
	print "Percent of library represented by final barcodes:", coverage

	return final_barcodes


def bootstrap_levenshtein(lib, n):
	"""
	This function calculates a reference Levenshtein distribution. It randomly
	picks two sequences from the reference sequences and calculates the distance
	to get a measure of how similar the library is.
	"""

	distances = []
	# bootstrap n times
	for i in range(0, n):
		# randomly grab two sequences with replacement
		string1 = random.choice(lib.keys())
		string2 = random.choice(lib.keys())

		distances.append(Levenshtein.distance(string1, string2))
	
	# take cutoff at 1% percentile
	cutoff = numpy.percentile(distances, 1)
	 
	return cutoff


if __name__ == '__main__':

	parser = argparse.ArgumentParser('Map barcodes to sequence')
	parser.add_argument('reads_file', help='sequence reads, accepted formats \
		FASTQ, FASTA or plain-text, one read per line')
	parser.add_argument('file_type', help='sequence read format, \
		specify fastq/fasta/txt')
	parser.add_argument('controls', help='Control library file if length shorter \
		than variant library, fasta format')
	parser.add_argument('controls_primer_len', type=int, help='primer length for control file')
	parser.add_argument('bc_loc', help='barcode location, specify start or end')
	parser.add_argument('bc_len', type=int, help='length of barcode')
	parser.add_argument('output_prefix', help='Name of output file prefix')
	parser.add_argument('--cutoff', type=int, help='user defined Levenshtein cutoff, \
		if not given then empirically determined (bootstrapped')

	args = parser.parse_args()

	reads_file = args.reads_file
	file_type = args.file_type
	bc_loc = args.bc_loc
	bc_len = args.bc_len
	output_prefix = args.output_prefix

	print ("Mapping controls...")
	if file_type == 'fastq':
		reads = fastq_reader(reads_file)
	elif file_type == 'fasta':
		# returns dictionary, just need reads (values) and not names (keys)
		reads = fasta_reader(reads_file).values()
	elif file_type == 'txt':
		reads = [line.strip() for line in open(reads_file)]
	else:
		raise Exception("Please specify either fastq, fasta or txt (raw one read per line)")

	print "Number of reads:", len(reads)
	
	# controls = library_reader(
	# 	filename=args.controls,
	# 	primer_len=args.controls_primer_len,
	# 	rev_complement=True,
	# 	format='tab')

	controls_fasta = fasta_reader(args.controls)
	# flip mapping so sequence is key and name is value, add reverse complement
	controls = {controls_fasta[x] : x for x in controls_fasta}
	controls.update({reverse_complement(controls_fasta[x]) : x for x in controls_fasta})

	print "Running bbmap..."
	# let's use bbmap to align the controls, much faster
	DEVNULL = open(os.devnull, 'w')
	command = shlex.split('bbmap.sh in={} ref={} out={} \
    	maxindel=200 nodisk=t noheader=t minid=0.20 threads=35'.format(
    		args.reads_file, args.controls, args.output_prefix + '_bbmap.sam'))

	bbmap = subprocess.Popen(command,
                         stdout=subprocess.PIPE,
                         stdin=subprocess.PIPE,
                         stderr=DEVNULL)

	bbmap.communicate()

    # parse SAM file
	aligned_reads, alignments = sam_parser(args.output_prefix + '_bbmap.sam')


	# grab barcodes that map to a aligned sequence
	if bc_loc == 'start':
		barcodes = [read[:bc_len] for read in aligned_reads]
	elif bc_loc == 'end':
		barcodes = [read[-bc_len:] for read in aligned_reads]


	print "Number of unique barcodes for aligned reads: ", len(set(barcodes))
	print "Filter by barcode frequency..."

	# Count barcodes 
	barcode_counts = dict(Counter(barcodes))

	# Throw out barcodes that appear 1 or 2 times, sequencing errors
	barcodes_clean = {x : barcode_counts[x] for x in barcode_counts if barcode_counts[x] > 2}
	print "Number of barcodes > 2:", len(barcodes_clean)

	variant_map, barcode_map = mapping(
		barcodes_clean, 
		aligned_reads, 
		reads,
		alignments,
		controls, 
		bc_loc, 
		bc_len)

	# bootstrap reference sequences to get a reference Levenshtein distribution 
	# to determine cutoff
	if args.cutoff:
		cutoff = args.cutoff
	else:
		print "Bootstrapping reference sequences to obtain cutoff...", 
		cutoff = bootstrap_levenshtein(controls, 10000)
		print "cutoff is Levenshtein distance ", cutoff

	print "Filtering and writing resultbs..."
	final_barcodes = filter_barcodes(
		barcode_map, 
		controls, 
		alignments,
		cutoff,
		bc_loc,
		bc_len,
		output_prefix + '_bc_map.txt')

	print "Number of final barcodes: ", len(final_barcodes)
