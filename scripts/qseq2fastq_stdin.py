'''
This script converts qseq files to FASTQ format. Takes as input
files piped in from stdin and writes output to stdout.
'''

import sys

if __name__ == '__main__':

	for line in sys.stdin:
		fields = line.strip().split('\t')
		if fields[10] == '1': # remove reads that don't pass filter which == 0
			header = '@' + ':'.join(fields[0:6]) + '#' + fields[6] + '/' + fields[7]
			print(header)
			# replace uncalled bases with N
			seq = fields[8].replace('.', 'N')
			print(seq) 
			print('+')
			print(fields[9])
