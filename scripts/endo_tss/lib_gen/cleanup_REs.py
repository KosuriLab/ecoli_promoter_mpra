import sys
import argparse
import itertools

def fasta_reader(fasta):
	'''
	Read in FASTA into dictionary format where keys are headers and 
	sequences are values
	'''
	promoters = {}

	while True:
		header = fasta.readline().strip()
		if not header: break
		seq = fasta.readline().strip()
		promoters[header] = seq

	print "Initial number of sequences:", len(promoters)
	
	return promoters


def remove_bad_promoters(promoters, res, standard, print_bad=False):
	good_promoters = {}
	bad_promoters = {}

	for x in promoters:
		seq = promoters[x]
		# trim primers
		trimmed = seq[24:-24]
		good = True
		for re in res:
			num_matches = trimmed.count(res[re])
			if re in standard:
				if num_matches > 1:
					good = False
					print x + " has RE site " + re + " more than once"
			else:
				if num_matches > 0:
					good = False
					print x + " has RE site " + re + " at least once"
		if good:
			good_promoters[x] = seq
		else:
			bad_promoters[x] = seq

	print "Removed", len(bad_promoters), "sequences."
	print "Final number of oligos:", len(good_promoters)

	return good_promoters	


def write_output(output_file, output_type, promoters):

	output = open(output_file, 'w')
	output_type = output_type

	if output_type == 'fasta':
		for x in good_promoters:
			output.write(x+'\n')
			output.write(good_promoters[x]+'\n')
	elif output_type == 'tab':
		for x in good_promoters:
			output.write(x+'\t')
			output.write(good_promoters[x]+'\n')
	else:
		print "Please specify output format. Valid arguments: tab or fasta"

	output.close()


if __name__ == '__main__':

	parser = argparse.ArgumentParser(description='Remove sequences that contain RE sites')
	parser.add_argument('promoter_filename', help='FASTA file of promoter sequences')
	parser.add_argument('re_file', help='FASTA file of restriction sites')
	parser.add_argument('output_file', help='Name of output file')
	parser.add_argument('output_type', help='Output format. tab or fasta')
	args = parser.parse_args()

	promoter_filename = args.promoter_filename
	promoter_file = open(promoter_filename, 'r')
	promoters = fasta_reader(promoter_file)

	re_file = open(args.re_file, 'r')
	REs = {}
	while True:
		lines = map(str.strip, list(itertools.islice(re_file, 2)))
		if not lines:
			break
		# strip '>' from name
		REs[lines[0][1:]] = lines[1]

	
	good_promoters = remove_bad_promoters(promoters, REs, standard=['XhoI', 'NheI'])

	write_output(args.output_file, args.output_type, good_promoters)
