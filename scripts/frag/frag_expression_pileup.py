import argparse
import numpy as np

def in_range(x, start, end):
	if x >= start and x <= end:
		return True
	else:
		return False


def pileup(frags, start_position, end_position, num_no_overlap, outfile_name):
	# frags must be sorted
	outfile = open(outfile_name, 'w')
	# write wig header
	outfile.write('variableStep chrom=U00096.2\n')
	current_frags = []
	# frag_pileup = []
	for i in range(start_position, end_position):
		
		if i % 50000 == 0:
			print "Position ", i, "..."

		# don't remove fragment until current position past end of first fragment
		if (i + 1) > frags[0][2]:
			frags.pop(0)
		
		no_overlap = 0
		for frag in frags:
			start = frag[1]
			end = frag[2]
			# fragment coordinates are 1-based
			overlap = in_range(i + 1, start, end)
			if overlap:
				current_frags.append(frag)
			else:
				no_overlap += 1
			if no_overlap >= num_no_overlap:
				break


		if len(current_frags) == 0:
			mean_exp = 0

		else:
			# take average expression of all overlapping frags
			mean_exp = round(np.mean([frag[0] for frag in current_frags]), 2)

		current_frags = []

		# frag_pileup.append((i, mean_exp))
		outfile.write(str(i+1) + '\t' + str(mean_exp) + '\n')

		# reset list
		current_frags = []

		if len(frags) == 0:
			break

	outfile.close()



if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('infile', help='tab-separated file of fragment expression data. \
		File must be sorted by start position and have a header.\
		Start must always be less than end. \
		4th field is expression, 10th field is start, 11th is end, 12th is strand')
	parser.add_argument('outfile', help='name of output file, wig format')
	parser.add_argument('num_no_overlap', type=int, help='Number of non-overlapping fragments to search before quitting')

	args = parser.parse_args()
	prefix = args.outfile
	num_no_overlap = args.num_no_overlap

	plus_frags = []
	minus_frags = []
	with open(args.infile) as infile:
		# read through header
		infile.readline()
		for line in infile:
			fields = line.strip().split()
			expression = float(fields[3])
			start = int(fields[9])
			end = int(fields[10])
			strand = fields[11]
			if strand == '+':
				plus_frags.append((expression, start, end, strand))
			else:
				# switch start and end so start is always less than end
				minus_frags.append((expression, start, end, strand))

		# make sure they are sorted
		plus_frags = sorted(plus_frags, key=lambda x: x[1])
		minus_frags = sorted(minus_frags, key=lambda x: x[1])
		
		print "Plus strand pileup..."
		pileup(plus_frags, 1, 4639310, num_no_overlap, 'plus_' + prefix)
		print "Minus strand pileup..."
		pileup(minus_frags, 1, 4639310, num_no_overlap, 'minus_' + prefix)

