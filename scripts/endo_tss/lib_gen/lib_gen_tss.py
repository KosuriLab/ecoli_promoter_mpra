'''
This script takes as input a tab-delimited text file containing the TSS in the 
first column and the strand in the second column. The first line must be a 
header. Coordinates for the input TSS file and for the genome must be the same. 
Reads in forward and reverse primer files. Chooses the strand with the 
lowest adenine content for synthesis purposes. Output is in FASTA or csv format.
'''

import argparse
import random


def parse_fasta(fasta_file):
	'''
	Parse a FASTA file of the genome into a single string.
	'''
	fasta = open(fasta_file, 'r')
	# read through header
	header = fasta.readline()
	genome = [line.strip() for line in fasta.readlines()]
	genome = ''.join(genome)
	return genome


def parse_tss(tss_file):
	'''
	Parses a three column, tab-delimited text file. The first column is the TSS, 
	the second is the strand information and the third is the source. 
	Assumes first line is a header. Returns a dictionary of dictionaries 
	representation of TSS information.
	'''

	tss = open(tss_file, 'r')
	# read through header
	header = tss.readline()

	lines = [line.strip() for line in tss.readlines()]
	fields = [line.split('\t') for line in lines]

	# format as dictionary of dictionaries. Makes code more transparent when
	# you access values by their names instead of index
	tss = {}
	index = 1
	for field in fields:
		position = int(field[0])
		strand = field[1]
		source = field[2]
		tss[index] = { 'position': position, 'strand': strand, 'source': source}
		index += 1

	return tss


def parse_primer_files(fwd_filename, rev_filename):
	'''
	Parse individual forward and reverse primer files and pair them up
	'''
	fwd_file = open(fwd_filename, 'r')
	fwd_primers = [line.strip() for line in fwd_file.readlines()
			   if not line.startswith('>')]

	rev_file = open(rev_filename, 'r')
	rev_primers = [line.strip() for line in rev_file.readlines()
			   if not line.startswith('>')]

	paired_primers = zip(fwd_primers, rev_primers)

	return paired_primers


def reverse_complement(seq):
	'''
	Returns the reverse complement of a sequence
	'''

	# reverse the sequence
	rev = seq[::-1]
	rc = ''

	for nt in rev:
		if nt == 'A':
			rc += 'T'
		elif nt == 'T':
			rc += 'A'
		elif nt == 'C':
			rc += 'G'
		else: # nt  == 'G'
			rc += 'C'
	
	return rc


def best_A_content(oligo):
	'''
	Choose the strand with the lowest A content because A's are harder to
	synthesize.
	'''
	rc_oligo = reverse_complement(oligo)

	oligo_As = sum( [1 for nt in oligo if nt == 'A'] )
	rc_As = sum( [1 for nt in rc_oligo if nt == 'A'] )

	if oligo_As < rc_As:
		final_oligo = oligo
	else:
		final_oligo = rc_oligo

	return final_oligo


def get_sequence(tss, strand, genome):
	'''
	Extracts sequence from the genome based on TSS and strand orientation. The
	sequence is from -120 to +30 relative to the TSS, giving 150 bp oligos. 
	The strand with the lowest A content is returned.
	'''

	if strand == '+':

		oligo_start = tss - 120
		if oligo_start < 0:
			oligo_start = 0

		# add +30 of UTR
		oligo_stop = oligo_start + 150
		if oligo_stop > len(genome):
			oligo_stop = len(genome) - 1
			oligo_start = oligo_stop - 150

		oligo = genome[oligo_start:oligo_stop]
	
	if strand == '-':
		oligo_start = tss - 30
		if oligo_start < 0:
			oligo_start = 0

		oligo_stop = oligo_start + 150
		if oligo_stop > len(genome):
			oligo_stop = len(genome) - 1
			oligo_start = oligo_stop + 150
		# get reverse complement
		oligo = reverse_complement(genome[oligo_start:oligo_stop])

	return oligo


def add_primers(primers, oligos):

	fwd_primer, rev_primer = primers.pop()
	oligo_with_primers = { x : fwd_primer + oligos[x] + rev_primer 
						   for x in oligos}

	# get strand of oligo with lowest A content
	best_oligos = { x : best_A_content(oligo_with_primers[x])
					for x in oligo_with_primers}

	return best_oligos


def seq_writer(seqs, tss, output_name, output_type):

	outfile = open(output_name, 'w')

	for x in seqs:

		position = str(tss[x]['position'])
		strand = tss[x]['strand']
		source = tss[x]['source']
		sequence = seqs[x]

		if output_type == 'fasta':
			outfile.write(','.join(['>TSS_' + str(x) + '_' + source, position, strand])+'\n')
			outfile.write(sequence+'\n')
		if output_type == 'tab':
			outfile.write(','.join(['>TSS_' + str(x) + '_' + source, position, strand])+'\t')
			outfile.write(sequence+'\n')

	outfile.close()


if __name__ == '__main__':

	parser = argparse.ArgumentParser(description='construct 150bp oligos, -120 to + 30, based on TSS')
	parser.add_argument('genome', help='FASTA file of genome')
	parser.add_argument('tss_list', help='Tab-delimited text file with header. One column for position, one for strand')
	parser.add_argument('fwd_primer', help='FASTA format of forward primers. Must be in same order as reverse file.')
	parser.add_argument('rev_primer', help='FASTA format of reverse primers. Must be in same order as forward file.')
	parser.add_argument('output_file', help='Name of output file.')
	parser.add_argument('output_type', help='Type of output. fasta or tab')

	args = parser.parse_args()

	# read in genome
	fasta_file = args.genome
	genome = parse_fasta(fasta_file)

	# read in tss info
	tss = parse_tss(args.tss_list)

	# store ID and sequence
	oligos = { x: get_sequence( tss[x]['position'], tss[x]['strand'], genome)
			   for x in tss}

	print "Library size:", len(oligos)

	# Read in primers
	primers = parse_primer_files(args.fwd_primer, args.rev_primer)

	# add primers and choose strand with best A content
	best_oligos = add_primers(primers, oligos)
	seq_writer(best_oligos, tss, args.output_file, args.output_type)



	


	

