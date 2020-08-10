"""
Generate bokeh plot for each region/operon specified in bed file
"""

import argparse
import ecoli_mpra_bokeh as ec_bokeh
from bokeh.io import export_svgs

def export_svg(p, output_folder, name):
	p.output_backend = 'svg'
	export_svgs(p, '{0}/{1}.svg'.format(output_folder, name))


if __name__ == '__main__':
	parser = argparse.ArgumentParser("programmatically generate bokeh .png plot\
		for each region specified in bed file")
	parser.add_argument('bed_file', help='filename of bed file, bed format')
	parser.add_argument('buffer', type=int, help="Number of bases to add on ends of region")
	parser.add_argument('output_folder', help='relative path to output folder')
	args = parser.parse_args()

	conditions = ['M9']
	nbuff = args.buffer
	output_folder = args.output_folder

	with open(args.bed_file) as infile:
		for line in infile:
			# initialize 
			fields = line.strip().split('\t')
			name = fields[3]
			strand = fields[-1]
			region_start = int(fields[1])
			region_end = int(fields[2])
			start = region_start - nbuff
			end = region_end + nbuff
			print(name+'...')
			# create plot elements
			src = ec_bokeh.make_pileup_dataset(conditions, start, end, ec_bokeh.colors)
			src_gene = ec_bokeh.make_region_genes(ec_bokeh.genes, start, end)
			src_tss = ec_bokeh.make_tss_arrow(ec_bokeh.endo_tss_lb, start, end, scaled=False)
			src_rna = ec_bokeh.make_rna_line(start, end, conditions)
			p = ec_bokeh.make_region_plot(src, src_gene, src_tss, src_rna, fixed_yaxis=3)
			export_svg(p, output_folder, name) 
