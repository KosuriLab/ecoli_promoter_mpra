'''
This script generates control sequences for the oligo library. The first set
of controls are the synthetic promoters from Sri's PNAS paper and the promoters
can be found in Dataset_S01 (http://www.pnas.org/content/110/34/14024?tab=ds).
The table was manually converted to CSV first. The sequences are shorter and a
constant stuffer sequence must be added between the primer site and another RE
site (XhoI in this case) so the complete 150bp sequence will be amplified but
the stuffer sequence will be cut. 

The second set of controls are regions of the genome that are not expected
to be under regulation. Choose sequences that are at least 200bp from a TSS (on 
either strand). It is acceptable if they lie within genes. Use one of these 
regions as a stuffer sequence for the PNAS promoters.

Author: Kimberly Insigne
Email:  kiminsigne@gmail.com
'''


import sys
import re
import argparse
import oligo_utils
import random
import subprocess


def extract_syn_sequences(filename):
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


def get_control_seqs(features, genome, n):
	'''
	Randomly find 150 bp sequences that the ends of the sequence are no closer 
	than 200bp from a TSS (on either strand). So if there is a gap of 550bp between
	TSS this is a good control sequence
	'''

	# get TSS sites
	tss = [features[x]['TSS'] for x in features if features[x]['TSS'] != '']
	# convert all from string to float
	tss = map(float, tss)
	# now convert to int
	tss = map(int, tss)
	tss.sort()

	good_regions = []
	for i in range(len(tss)-1):
		curr_pos = tss[i]
		next_pos = tss[i+1]

		if next_pos - curr_pos > 550:
			good_regions.append( (curr_pos, next_pos) )

	# randomly choose control sequences from these good regions
	control_regions = random.sample(good_regions, n)
	control_oligos = {}

	# get oligos from the center of the region
	for region in control_regions:
		start, stop = region
		center = int( start + ( (stop - start)/2 ) ) 
		oligo_start = center - 75
		oligo_stop = center + 75
		oligo = genome[oligo_start:oligo_stop]
		control_oligos[ 'neg_control_' + str(oligo_start) + ':' + str(oligo_stop)] = oligo

	return control_oligos


def add_stuffer(syn_promoters, control_seqs, fwd_primer, rev_primer):
	'''
	Add stuffer  + XhoI site + promoter to make the whole thing 148 bp
	'''
	xhoI = 'CTCGAG'
	seqs = list(control_seqs.values())
	stuffer = random.choice(seqs)

	stuffed_promoters = {}

	for x in syn_promoters:
		promoter = syn_promoters[x]
		to_add = 150 - len(promoter)
		new_promoter = fwd_primer[:20] + stuffer[:to_add-2] + xhoI + promoter + rev_primer
		stuffed_promoters[x] = new_promoter

	return stuffed_promoters


def best_A_content(oligo):
	'''
	Choose the strand with the lowest A content because A's are harder to
	synthesize.
	'''
	rc_oligo = oligo_utils.reverse_complement(oligo)

	oligo_As = sum( [1 for nt in oligo if nt == 'A'] )
	rc_As = sum( [1 for nt in rc_oligo if nt == 'A'] )

	if oligo_As < rc_As:
		final_oligo = oligo
	else:
		final_oligo = rc_oligo

	return final_oligo


def seq_writer(seqs, output_name, output_type):

	outfile = open(output_name, 'w')

	for x in seqs:

		if output_type == 'fasta':
			outfile.write('>' + x + '\n' + seqs[x] + '\n')
		if output_type == 'tab':
			outfile.write('>' + x + '\t' + seqs[x] + '\n')

	outfile.close()

		
if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Generate control sequences for library')
	parser.add_argument('syn_file', help='CSV of PNAS synthetic promoters')
	parser.add_argument('feature_file', help='TSV of features')
	parser.add_argument('genome', help='FASTA file of genome')
	parser.add_argument('output', help='filename of output file')

	args = parser.parse_args()

	syn_promoters = extract_syn_sequences(args.syn_file)
	features = oligo_utils.parse_features(args.feature_file)
	genome = oligo_utils.genome_reader(args.genome)

	control_seqs = get_control_seqs(features, genome, 500)

	# >skpp-5-F sk20mer-35738
	fwd_primer = 'TCCGACGGGGAGTATATACTCGAG'
	# >skpp-22-R sk20mer-39638
	rev_primer = 'GCTAGCGGCTTGATAGTTGCATTA'

	stuffed_promoters = add_stuffer(syn_promoters, control_seqs, fwd_primer, rev_primer)

	# add primers to control_seqs
	control_seqs_with_primers = { x : fwd_primer + control_seqs[x] + rev_primer 
							      for x in control_seqs }

	# combine sequences
	control_seqs_with_primers.update(stuffed_promoters)
	# choose strand with best A content
	best_oligos = { x : best_A_content(control_seqs_with_primers[x]) 
					for x in control_seqs_with_primers}

	print "Number of control sequences:", len(best_oligos)

	seq_writer(best_oligos, args.output, 'tab')
