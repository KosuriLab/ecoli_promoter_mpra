# This script calls "peaks" in our genome fragmentation data. A peak is defined
# as a continguous region with expression (normalized RNA/normalized DNA) above
# a user defined threshold. Input must be in wig format, separated by strand. 
# Regions within a user defined bp distance will be combined into a single region.

import argparse
from wiggelen import walk, fill
from copy import deepcopy
from collections import OrderedDict


def merge_regions(called_regions, merge_dist):
		merged_regions = []
		curr_regions = called_regions
		merged = True
		while merged:
			merged = False
			for i in range(0, len(curr_regions) - 1, 2): # step through two at a time
				curr_start, curr_end, curr_value = curr_regions[i]
				next_start, next_end, next_value = curr_regions[i + 1]
				if next_start - curr_end <= merge_dist:
					merged = True
					new_region = (curr_start, next_end, curr_value + next_value)
					merged_regions.append(new_region)
				else:
					merged_regions.append(curr_regions[i])
					merged_regions.append(curr_regions[i + 1])
			
			# add last region if odd
			if len(curr_regions) % 2 == 1:
				merged_regions.append(curr_regions[-1])

			curr_regions = deepcopy(merged_regions)
			merged_regions = []

		return curr_regions

# wig = open('../../processed_data/frag_peak_calling/plus_frag_pileup.wig')
def find_max(regions, wig):
	# find max within a given region. Assumes regions are ordered and non-overlapping
	# find max by walking through wig file
	# Let's also calculate a new score sum just to be extra sure it's right. The
	# old score sum was done during merging regions, so it could have been messed up
	regions_with_max = []

	start, end, score_sum = regions.pop(0)
	max_value = -1
	max_position = None
	new_score_sum = 0
	for chrom, position, value in fill(walk(wig)):
		if start <= position and position <= end:
			new_score_sum += value
			if value >= max_value:
				max_value = value
				max_position = position
		if position >= end and max_value != -1:
			# position is past region and max_value has been recorded
			region_with_max = (start, end, new_score_sum, max_value, max_position)
			regions_with_max.append(region_with_max)
			# grab new region and reset max
			if len(regions) == 0: # break when there are no more regions
				break
			start, end, score_sum = regions.pop(0)
			max_value = -1
			max_position = None

	return regions_with_max


def write_bed(called_regions, strand, outfile):
	chrom = 'U00096.2'
	with open(outfile, 'w') as output:
		for x in called_regions:
			x_str = map(str, x)
			# format: chrom, start, end, name, score
			name = '_'.join([x_str[0], x_str[1], strand])
			output.write('\t'.join([chrom, x_str[0], x_str[1], name, x_str[2], strand, x_str[3], x_str[4]]))
			output.write('\n')


def main(infile, threshold, merge_dist, min_width, strand, outfile):
	# list of tuples, (start, end, avg_exp)
	called_regions = []
	
	start = None
	end = None
	total_exp = 0

	# fill function steps through every position, returns None if position not 
	# in original wig file
	print("Calling preliminary regions...")
	wig = open(infile)
	for region, position, value in fill(walk(wig)):
		if start is None:
			# initialize start of new region to current position
			start = position
		if value is None:
			if total_exp > 0: # if a region already exists, end it
				called_regions.append((start, end, total_exp))
			# reset start, end, and total_exp
			start = None
			end = None
			total_exp = 0
		elif value < threshold: 
			if total_exp > 0: # if a region already exists, end it
				called_regions.append((start, end, total_exp))
			# reset start, end, and total_exp
			start = None
			end = None
			total_exp = 0
		elif value >= threshold: # value exceeds threshold, continue region
			total_exp += value
			end = position
	
	wig.close()	

	if total_exp != 0: # finished iterating but one last region
		called_regions.append((start, end, total_exp))

	print("Filtering out regions smaller than minimum width...")
	# filter out regions that are below minimum width
	filtered_regions = [x for x in called_regions if x[1] - x[0] + 1 >= min_width]
	
	print("Merging regions...")
	merged_regions = merge_regions(filtered_regions, merge_dist)

	# find max region and re-do score sum
	# open wig file again
	print("Finding region max and calculating total score...")
	regions_with_max = find_max(merged_regions, open(infile))

	write_bed(regions_with_max, strand, outfile)


if __name__ == '__main__':

	parser = argparse.ArgumentParser()
	parser.add_argument('infile', help='strand specific wig file')
	parser.add_argument('threshold', type=float, help="expression threshold for peak calling")
	parser.add_argument('merge_dist', type=int, help='merge regions within n bp')
	parser.add_argument('min_width', type=int, help='minimum width of peak')
	parser.add_argument('strand', help='strand')
	parser.add_argument('outfile', help='name of output file, bed format')

	args = parser.parse_args()
	infile = args.infile
	threshold = args.threshold
	merge_dist = args.merge_dist
	min_width = args.min_width
	strand = args.strand
	outfile = args.outfile

	main(infile, threshold, merge_dist, min_width, strand, outfile)

	



