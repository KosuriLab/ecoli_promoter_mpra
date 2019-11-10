import argparse
import os
import call_peaks
import numpy as np


def count_overlap(tss_file, peak_file):
	# by default, if an overlap is found, reports shared interval between two
	# overlapping features
	# write original entry in B for each overlap
	cmd = 'bedtools intersect -wb -a ' + tss_file + ' -b ' + peak_file + ' > overlap.bed'
	os.system(cmd)
	# get unique peak overlaps (unique B entries)
	cmd = 'cut -f 7-12 overlap.bed | sort | uniq > overlap_unique.bed'
	os.system(cmd)
	num_overlap = sum(1 for line in open('overlap_unique.bed'))
	return num_overlap


def count_overlap_negative(neg_file, peak_file, pos_file):
	# for negatives, do not count as negative if it overlaps with at least one 
	# positive
	cmd = 'bedtools intersect -wb -a ' + neg_file + ' -b ' + peak_file + ' > overlap_negative.bed'
	os.system(cmd)
	# get unique peak overlaps (unique B entries)
	cmd = 'cut -f 7-12 overlap_negative.bed | sort | uniq > overlap_negative_unique.bed'
	os.system(cmd)
	# overlap with positive, get unique
	cmd = 'bedtools intersect -wb -a ' + pos_file + ' -b overlap_negative_unique.bed | cut -f 7-12 | sort | uniq > overlap_negative_unique_pos.bed'
	os.system(cmd)
	# do not count peaks that overlap with positives as negatives
	num_overlap = sum(1 for line in open('overlap_negative_unique.bed')) - sum(1 for line in open('overlap_negative_unique_pos.bed'))
	return num_overlap


if __name__ == '__main__':

	parser = argparse.ArgumentParser()
	parser.add_argument('pos_file', help='bed file of positive TSSs')
	parser.add_argument('neg_file', help='bed file of negative TSSs')
	parser.add_argument('plus_wig', help='wig file for plus strand')
	parser.add_argument('minus_wig', help='wig file for minus strand')
	parser.add_argument('min_threshold', type=float, help='minimum threshold')
	parser.add_argument('max_threshold', type=float, help='maximum threshold')
	parser.add_argument('threshold_step', type=float, help='step size for threshold')
	parser.add_argument('merge_dist', type=int, help='distance between adjacent peaks to merge')
	parser.add_argument('min_width', type=int, help='min width of peak')

	args = parser.parse_args()
	pos_file = args.pos_file
	neg_file = args.neg_file
	min_threshold = args.min_threshold
	max_threshold = args.max_threshold
	step = args.threshold_step
	merge_dist = args.merge_dist
	min_width = args.min_width

	# separate positive and negative files by strand
	cmd = 'awk \'{if ($6 == \"+\") print $0}\' ' + pos_file + ' > pos_file_plus.bed'
	os.system(cmd)

	cmd = 'awk \'{if ($6 == \"-\") print $0}\' ' + pos_file + ' > pos_file_minus.bed'
	os.system(cmd)

	cmd = 'awk \'{if ($6 == \"+\") print $0}\' ' + neg_file + ' > neg_file_plus.bed'
	os.system(cmd)

	cmd = 'awk \'{if ($6 == \"-\") print $0}\' ' + neg_file + ' > neg_file_minus.bed'
	os.system(cmd)


	outfile = open('optimize_peak_call_results.txt', 'w')
	# write header
	outfile.write('\t'.join(['threshold', 'merge_dist', 'min_width', 
		'num_peaks_tss_overlap_plus', 'num_peaks_tss_overlap_minus',
		'num_peaks_plus', 'num_peaks_minus',
		'plus_positive_overlap', 'plus_negative_overlap',
		'minus_positive_overlap', 'minus_negative_overlap', ]) + '\n')

	for threshold in np.arange(min_threshold, max_threshold, step):
		print "Threshold: ", threshold
		# plus strand
		call_peaks.main(args.plus_wig, threshold, merge_dist, min_width, '+', 'tmp_plus_peaks.bed')
		# minus strand
		call_peaks.main(args.minus_wig, threshold, merge_dist, min_width, '-', 'tmp_minus_peaks.bed')

		num_peaks_plus = sum(1 for line in open('tmp_plus_peaks.bed'))
		num_peaks_minus = sum(1 for line in open('tmp_minus_peaks.bed'))

		# calculate overlap between peak calls and TSSs
		plus_positive_overlap = count_overlap('pos_file_plus.bed', 'tmp_plus_peaks.bed')
		plus_negative_overlap = count_overlap_negative('neg_file_plus.bed', 'tmp_plus_peaks.bed', 'pos_file_plus.bed')

		minus_positive_overlap = count_overlap('pos_file_minus.bed', 'tmp_minus_peaks.bed')
		minus_negative_overlap = count_overlap_negative('neg_file_minus.bed', 'tmp_minus_peaks.bed', 'pos_file_minus.bed')

		num_peaks_tss_overlap_plus = plus_positive_overlap + plus_negative_overlap
		num_peaks_tss_overlap_minus = minus_positive_overlap + minus_negative_overlap

		info = [threshold, merge_dist, min_width,
		num_peaks_tss_overlap_plus, num_peaks_tss_overlap_minus,
		num_peaks_plus, num_peaks_minus,
		plus_positive_overlap, plus_negative_overlap, minus_positive_overlap, minus_negative_overlap]
		info = map(str, info)

		outfile.write('\t'.join(info) + '\n')


	outfile.close()
	# cmd = 'rm -f tmp_plus_peaks.bed tmp_minus_peaks.bed pos_file_plus.bed pos_file_minus.bed neg_file_plus.bed neg_file_minus.bed overlap*.bed'
	# os.system(cmd)