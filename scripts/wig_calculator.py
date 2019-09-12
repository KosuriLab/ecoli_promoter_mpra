import argparse
from wiggelen import walk, fill


def add(x, y):
	return x + y


def subtract(x, y):
	return x - y


def division(x, y):
	if y == 0:
		return 0
	else:
		return float(x) / float(y)


if __name__ == '__main__':

	parser = argparse.ArgumentParser()
	parser.add_argument('file1', help='name of first wig file')
	parser.add_argument('file2', help='name of second wig file')
	parser.add_argument('operation', help='''Mathematical operation. add, subtract
		or divide''')
	parser.add_argument('output', help='Name of output file')
	parser.add_argument('--median_norm', action='store_true')

	args = parser.parse_args()

	file1 = args.file1
	file2 = args.file2
	operation = args.operation

	print "Reading file 1..."
	wig1_dict = {position : value for region, position, value in fill(walk(open(file1)))}
	print "Reading file 2.."
	wig2_dict = {position : value for region, position, value in fill(walk(open(file2)))}

	min_position = min(min(wig1_dict.keys()), min(wig2_dict.keys()))
	max_position = max(max(wig1_dict.keys()), max(wig2_dict.keys()))

	print "Calculating..."
	if operation == 'add':
		output_dict = {i : add(wig1_dict.get(i, 0), wig2_dict.get(i, 0)) for i in range(min_position, max_position+1)}
	elif operation == 'subtract':
		output_dict = {i : subtract(wig1_dict.get(i, 0), wig2_dict.get(i, 0)) for i in range(min_position, max_position+1)}
	elif operation == 'division':
		output_dict = {i : division(wig1_dict.get(i, 0), wig2_dict.get(i, 0)) for i in range(min_position, max_position+1)}
	else:
		raise ValueError('invalid operation')

	outfile = open(args.output, 'w')
	outfile.write("variableStep chrom=U00096.2\n")

	for i in range(min_position, max_position+1):
		outfile.write(str(i) + '\t' + str(output_dict[i]) + '\n')


