import sys
import argparse

parser = argparse.ArgumentParser()
parser.add_argument(infile, help='input file')
parser.add_argument(outfile, help='output file')
args = parser.parse_args()

infile = open(args.infile)
outfile = open(args.outfile, 'w')

for line in infile.readlines():
	fields = line.strip().split('\t')
	header = fields[0]
	seq = fields[1]

	outfile.write(header+'\n')
	outfile.write(seq+'\n')

outfile.close()