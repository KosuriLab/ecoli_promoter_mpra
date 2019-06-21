import argparse
import pandas as pd
from wiggelen import walk, fill


def wig_lookup(start, end, wig):

	wig_sum = sum([wig.get(i, 0) for i in range(start, end)])
	return wig_sum


if __name__ == '__main__':

	parser = argparse.ArgumentParser()
	parser.add_argument('bed', help='bed file')
	parser.add_argument('plus_wig', help='wig file for plus strand')
	parser.add_argument('minus_wig', help='wig file for minus strand')
	parser.add_argument('output_name', help='name of output file')

	args = parser.parse_args()

	plus_wig_dict = {position : value for region, position, value in fill(walk(open(args.plus_wig)))}
	minus_wig_dict = {position : value for region, position, value in fill(walk(open(args.minus_wig)))}

	outfile = open(args.output_name, 'w')

	with open(args.bed) as infile:
		for line in infile:
			fields = line.strip().split()
			chrom = fields[0]
			start = int(fields[1])				
			end = int(fields[2])
			name = fields[3]
			strand = fields[5]

			if strand == '+':
				wig_sum = wig_lookup(start, end, plus_wig_dict)
			elif strand == '-':
				wig_sum = wig_lookup(start, end, minus_wig_dict)
			else:
				wig_sum = -1

			outfile.write('\t'.join(map(str, [chrom, start, end, name, wig_sum, strand])) + '\n')

	outfile.close()

