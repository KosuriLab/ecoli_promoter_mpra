# This script parses the GSE52059 transcriptome annotation file
import sys

annot = open(sys.argv[1])
output = open(sys.argv[2], 'w')

# read through header
annot.readline()

while True:

	line = annot.readline()
	if not line: 
		break

	fields = line.strip().split()
	if len(fields) == 2: # if this is not a comment line, e.g /promoter='primary'
		name = fields[0]
		location = fields[1]
	else:
		continue

	if name == 'promoter':
		# parse the location field
		if location.startswith('complement'):
			strand = '-'
			# format: complement(100..100)
			# split by '(' then by '..'
			position = location.split('(')[1].split('..')[0]
		else:
			strand = '+'
			# format: 100..100 , split by '..'
			position = location.split('..')[0]
		# get promoter type information on next line
		type_info = annot.readline().strip().split()[0]
		# format: /promoter="primary" , split by "
		promoter_type = type_info.split('\"')[1]

		# four possible promoter types: primary, secondary, internal and antisense
		# encoded as 1 if it is that type, 0 otherwise. This matches the notation
		# of the other Storz TSS file. If promoter is primary, type will be
		# encoded as 1,0,0,0. 
		types = ['0', '0', '0', '0']

		if promoter_type == 'primary':
			types[0] = '1'
		elif promoter_type == 'secondary':
			types[1] = '1'
		elif promoter_type == 'internal':
			types[2] = '1'
		elif promoter_type == 'anti-sense':
			types[3] = '1'

		types = ','.join(types)

		output.write('\t'.join([position, strand, types]) + '\n')


