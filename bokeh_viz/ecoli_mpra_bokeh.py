from bokeh.plotting import figure, output_notebook, show, ColumnDataSource, output_file
from bokeh.models.plots import Plot
from bokeh.models import Legend, LegendItem
from bokeh.models.tools import HoverTool, BoxSelectTool, BoxZoomTool, PanTool, WheelZoomTool, SaveTool, ResetTool
from bokeh.models.tickers import FixedTicker
from bokeh.models.widgets import CheckboxGroup, RangeSlider, Tabs, TextInput, Button, RadioGroup, RadioButtonGroup, Select, Paragraph, Div
from bokeh.layouts import column, row, WidgetBox, Spacer, gridplot
from bokeh.models import Panel, Range1d, LinearAxis
from bokeh.io import show, curdoc, export_svgs
import os
import pandas as pd
import numpy as np
from math import pi
from copy import deepcopy
import sys

# change to script directory
abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
os.chdir(dname)


global frag_plus_pileup
global frag_minus_pileup
global colors
global genes
global rna_plus
global rna_minus

NEW_MIN = 1
NEW_MAX = 1000


def norm_pileup(pileup):
    # exclude data below median expression. However, for graphing we need a continuous line
    # for all positions. We will set all data below median to the median value, and then scale
    # data so minimum value is 1. On log-scale, log(1) = 0 so any data <= median will be set to
    # zero on plot
    median = np.median(pileup.expression)
    norm_pileup = deepcopy(pileup)
    norm_pileup.expression = np.where(pileup.expression <= median,
                                     median, pileup.expression)
    scalar = 1.0/median
    norm_pileup.expression = norm_pileup.expression * scalar
    # make sure everything is equal to 1, due to decimal multiplication some may be 0.9999
    norm_pileup.expression = np.where(norm_pileup.expression < 1,
                                     1, norm_pileup.expression)
    return norm_pileup


def min_max_scale(expression, new_min, new_max, max_value=None, median_norm=None):

    if max_value:
        # cap value at max so the scaling won't be skewed
        expression = np.where(expression > max_value, max_value, expression)

    if median_norm:
        # subtract median
        expression = expression - median_norm
        # if negative, set to 0
        expression = np.where(expression < 0, 0, expression)

    # first min-max normalize to [0, 1]
    min_exp = min(expression)
    max_exp = max(expression)

    norm = (expression - min_exp) / float(max_exp - min_exp)
    norm_scaled = (norm * (new_max - new_min)) + new_min

    return norm_scaled


def make_pileup_dataset(conditions, start, end, colors, log_transform=True):
    
    # subset relevant condition, conditions must be list
    plus_pileup_subset = frag_plus_pileup_norm[(frag_plus_pileup_norm.condition.isin(conditions)) &
                                         (frag_plus_pileup_norm.position.isin(range(start,end)))]
    
    minus_pileup_subset = frag_minus_pileup_norm[(frag_minus_pileup_norm.condition.isin(conditions)) &
                                           (frag_minus_pileup_norm.position.isin(range(start,end)))]
    
    plus_pileup_subset = plus_pileup_subset.merge(colors[colors.strand == '+'], on='condition', how='left')
    minus_pileup_subset = minus_pileup_subset.merge(colors[colors.strand == '-'], on='condition', how='left')
    
    # patch accepts list of lists, make one list for each combination of strand and condition
    # x-axis position will be same for every condition and strand, define coordinates for edges at the top (pileup line) and the 
    position_unique = plus_pileup_subset.sort_values('position').position.unique()
    position = np.hstack((position_unique, position_unique[::-1]))
    # bottom edge along zero (positions in reverse to create closed shape)
    bottom_edge = np.zeros(len(position))
    # take negative of reverse coverage to display on graph
    counts = []
    position_all = []
    condition_colors = []
    for condition in conditions:
        # plus strand
        position_all.append(position)
        if log_transform:
            counts.append(np.hstack((np.log10(plus_pileup_subset.expression[plus_pileup_subset.condition == condition]),
                               bottom_edge)))
        else:
            counts.append(np.hstack((plus_pileup_subset.expression[plus_pileup_subset.condition == condition],
                                   bottom_edge)))
        condition_colors.append(colors.color[(colors.condition == condition) & (colors.strand == '+')].to_string(index=False))
        
        # minus strand
        position_all.append(position)
        if log_transform:
            counts.append(np.hstack((-np.log10(minus_pileup_subset.expression[minus_pileup_subset.condition == condition]),
                               bottom_edge)))
        else:
            counts.append(np.hstack((-minus_pileup_subset.expression[minus_pileup_subset.condition == condition],
                                   bottom_edge)))
        condition_colors.append(colors.color[(colors.condition == condition) & (colors.strand == '-')].to_string(index=False))
    
    patches_src = pd.DataFrame({'position' : position_all, 
                               'count' : counts,
                               'color' : condition_colors})    
    
    return ColumnDataSource(patches_src)


def make_region_genes(genes, start, end):
    # grab genes within coordinates, with buffer of 100bp
    start_buffer = start - 100
    end_buffer = end + 100
    region_genes = genes[(genes.start >= start_buffer) & (genes.end <= end_buffer)]

    # center is midpoint of gene
    gene_center = region_genes.start + abs(region_genes.end - region_genes.start)/2.0
    gene_center = gene_center.tolist()
    gene_color = '#e2b306'
    gene_width = region_genes.end - region_genes.start
    gene_width = gene_width.tolist()

    src_gene = ColumnDataSource(data=dict(
    gene_center = gene_center,
    # y-center is on x-axis, i.e y = 0
    gene_center_y = np.zeros(len(region_genes)),
    gene_color = [gene_color] * len(region_genes),
    gene_width = gene_width,
    # set center of triangle to start or end depending on strand
    tri_x = [region_genes.end.iloc[i] if region_genes.strand.iloc[i] == '+' else region_genes.start.iloc[i] \
           for i in range(len(region_genes))],
    angle = [-90 if strand == '+' else 90 for strand in region_genes.strand],
    gene_name = region_genes.name.tolist()))
    
    return src_gene


def make_tss_arrow(tss, start, end, threshold=1, active_color='#2e6eb7', inactive_color='grey', width=8, arrow_length=50):

    # grab TSS within start and end
    tss_region = tss[tss.tss_position.isin(range(start, end))]
    
    # list of lists. Need four quad coordinates: top, bottom, left, right
    # top and bottom are y-coordinates, left and right are x-coordinates
    
    # draw vertical rectangle to indicate TSS position. Above axis indicates positive strand TSS
    # and below axis indicates negative strand. We'll keep it as a list of lists [top, bottom, left, right]
    # and separate into different lists in ColumnDataSource
    quad_coords = []
    
    # draw horizontal segment to connect to arrow, [x0, y0, x1, y1]
    seg_coords = []
    
    # center of triangle is segment endpoint. Set direction of triangle based on strand
    arrow_angle = []
    
    # assign color based on TSS expression above or below threshold
    tss_color = []
    
    for i in range(len(tss_region)):
        
        expn = tss_region.RNA_exp_sum_ave.iloc[i]
        position = tss_region.tss_position.iloc[i]
        
        if expn < threshold:
            tss_color.append(inactive_color)   
        else:
            tss_color.append(active_color)
        
        if tss_region.strand.iloc[i] == '+':
            # top, bottom, left, right
            quad_coords.append([expn, 0, position, position + width])
            # draw line to the right: x0, y0, x1, y1
            seg_coords.append([position, expn, position + arrow_length, expn])
            arrow_angle.append(-90)
        elif tss_region.strand.iloc[i] == '-':
            quad_coords.append([0, -expn, position, position + width])
            # draw line to the left
            seg_coords.append([position, -expn, position - arrow_length, -expn])
            arrow_angle.append(90)
        else:
            raise Exception('Invalid strand')
    
    src_tss = ColumnDataSource(data=dict(
    top = [quad_coords[i][0] for i in range(len(quad_coords))],
    bottom = [quad_coords[i][1] for i in range(len(quad_coords))],
    left = [quad_coords[i][2] for i in range(len(quad_coords))],
    right = [quad_coords[i][3] for i in range(len(quad_coords))],
    x0 = [seg_coords[i][0] for i in range(len(seg_coords))],
    y0 = [seg_coords[i][1] for i in range(len(seg_coords))],
    x1 = [seg_coords[i][2] for i in range(len(seg_coords))],
    y1 = [seg_coords[i][3] for i in range(len(seg_coords))],
    angle = arrow_angle,
    color = tss_color))
    
    return src_tss


def make_rna_line(start, end, conditions):

    # subset for relevant condition and position
    rna_plus_subset = rna_plus[(rna_plus.condition.isin(conditions)) &
                                           (rna_plus.position.isin(range(start, end)))]
    rna_minus_subset = rna_minus[(rna_minus.condition.isin(conditions)) &
                                           (rna_minus.position.isin(range(start, end)))]
    rna_minus_subset.head()

    # convert expression to log and flip minus expression to negative
    src_rna = ColumnDataSource(data=dict(
    x_plus = rna_plus_subset.position.tolist(),
    y_plus = np.log10(rna_plus_subset.expression.tolist()),
    x_minus = rna_minus_subset.position.tolist(),
    y_minus = -(np.log10(rna_minus_subset.expression.tolist()))))

    return src_rna


def make_region_plot(src, src_gene, src_tss, src_rna):
    '''
    Construct pileup plot based on src
    '''

    # output_file(html_output)

    # flatten list of lists in count column of src, find max value of absolute expression
    frag_count_range = max(map(abs, [x for count in src.data['count'] for x in count]))
    if len(src_rna.data['y_plus']) > 0:
        rna_count_range = max(max(src_rna.data['y_plus']), max(abs(src_rna.data['y_minus'])))
    else: 
        rna_count_range = 0

    count_range = max(frag_count_range, rna_count_range)
    
    # draw blank figure of correct size with tools
    p = figure(y_range=(-(count_range+1), count_range+1), plot_width=900, plot_height=700, 
               tools=[BoxSelectTool(), BoxZoomTool(), PanTool(), WheelZoomTool(), 
                      SaveTool(), ResetTool()],
               toolbar_location='above')

    legends = []
    
    # format axis and colors
    p.xaxis.axis_label = 'position'
    p.xaxis.major_label_orientation = pi/4
    p.xaxis[0].formatter.use_scientific = False
    # p.xaxis[0].ticker=FixedTicker(ticks=range(start, end, 100))
    p.yaxis.axis_label = 'log normalized activity'
        
    patches = p.patches(source=src, xs='position', ys='count', 
        fill_color='color', line_color=None, alpha=0.50)
    legends.append(LegendItem(label='promoter activity (plus strand)', 
        renderers=[patches], index=0))
    legends.append(LegendItem(label='promoter activity (minus strand)', 
        renderers=[patches], index=1))
    

    # draw RNA lines
    if len(src_rna.data['x_minus']) > 0:
        plus_line = p.line(x='x_plus', y='y_plus', line_color='#528ecb', 
            line_width=2, source=src_rna)
        legends.append(LegendItem(label='RNA-seq (plus strand)',
            renderers=[plus_line], index=0))
        minus_line = p.line(x='x_minus', y='y_minus', line_color='#ef8137',
            line_width=2, source=src_rna)
        legends.append(LegendItem(label='RNA-seq (minus strand)',
            renderers=[minus_line], index=1))


    # add second y-axis for TSS strength
    max_tss = max(map(abs, src_tss.data['y0']))
    p.extra_y_ranges = {'tss' : Range1d(start=-(max_tss * 4), end=(max_tss * 4))}
    p.add_layout(LinearAxis(y_range_name='tss', axis_label='TSS expression'), 'right')
    
    # draw vertical rectangle
    p.quad(top='top', bottom='bottom', left='left', right='right', color='color',
           source=src_tss, y_range_name='tss')
    # draw horizontal line for arrow
    p.segment(x0='x0', y0='y0', x1='x1', y1='y1', color='color',
              source=src_tss, line_width=4,  y_range_name='tss')
    # center of triangle is endpoint of segment
    tri = p.triangle(x='x1', y='y1', size=9, angle='angle', angle_units='deg',
        color='color', source=src_tss, y_range_name='tss')
    legends.append(LegendItem(label='inactive TSS', renderers=[tri], index=9))
    legends.append(LegendItem(label='active TSS', renderers=[tri], index=0))
    

    # plot genes
    p.rect(x='gene_center', y='gene_center_y', width='gene_width', color='gene_color', 
        height=10, height_units='screen', source=src_gene)
    p.triangle(x='tri_x', y=0, size=20, angle='angle', angle_units='deg',
                     fill_color='gene_color', line_color=None, source=src_gene)
    p.text(x='gene_center', y='gene_center_y', text='gene_name', text_color='black',
          text_align='center', text_baseline='middle', text_font_size='10pt', source=src_gene)
    

    p.add_layout(Legend(items=legends), 'below')

    return p
        

def update_region_plot(attr, old, new):
    '''
    Update plot based on interactive widgets
    '''
    # get list of conditions for graph
    conditions_to_plot = [condition_selection.labels[i] for i in condition_selection.active]
    

    position_start = int(start.value)
    position_end = int(end.value)
    
    # make new subset based on selected conditions
    new_src = make_pileup_dataset(conditions_to_plot, position_start, position_end, colors)
    new_src_gene = make_region_genes(genes, position_start, position_end)
    new_src_tss = make_tss_arrow(endo_tss_lb, position_start, position_end)
    new_src_rna = make_rna_line(position_start, position_end, conditions_to_plot)
    
    # update source used in patch glyphs
    src.data.update(new_src.data)
    # update source used in gene rectangle and triangle glyphs
    src_gene.data.update(new_src_gene.data)
    # update source for TSS arrows
    src_tss.data.update(new_src_tss.data)
    # update source for RNA
    src_rna.data.update(new_src_rna.data)


def toggle_genes(attr, old, new):

    if gene_button.value == 'show': # show genes
        # grab genes at current input position, update
        new_src_gene = make_region_genes(genes, int(start.value), int(end.value))
        src_gene.data.update(new_src_gene.data)
    if gene_button.value == 'hide': # hide genes
        # no genes
        new_src_gene = make_region_genes(genes, 0, 0)
        src_gene.data.update(new_src_gene.data)


def toggle_tss(attr, old, new):

    if tss_button.value == 'show':
        new_src_tss = make_tss_arrow(endo_tss_lb, int(start.value), int(end.value))
        src_tss.data.update(new_src_tss.data)
    if tss_button.value == 'hide':
        new_src_tss = make_tss_arrow(endo_tss_lb, 0, 0)
        src_tss.data.update(new_src_tss.data)


def toggle_rna(attr, old, new):

    if rna_button.value == 'show':
        new_src_rna = make_rna_line(int(start.value), int(end.value), 
            [condition_selection.labels[i] for i in condition_selection.active])
        src_rna.data.update(new_src_rna.data)
    if rna_button.value == 'hide':
        new_src_rna = make_rna_line(0, 0, [condition_selection.labels[i] for i in condition_selection.active])
        src_rna.data.update(new_src_rna.data)



def update_axis():
    # flatten list of lists in count column of src, find max value of absolute expression
    frag_count_range = max(map(abs, [x for count in src.data['count'] for x in count]))
    if len(src_rna.data['y_plus']) > 0:
        rna_count_range = max(max(src_rna.data['y_plus']), max(abs(src_rna.data['y_minus'])))
    else: 
        rna_count_range = 0

    count_range = max(frag_count_range, rna_count_range)

    p.y_range.start = -(count_range + 1)
    p.y_range.end = count_range + 1

    # update second y-axis
    max_tss = max(map(abs, src_tss.data['y0']))
    p.extra_y_ranges['tss'].start = max_tss * -4
    p.extra_y_ranges['tss'].end = max_tss * 4


def add_bp(n=500):

    start.value = str(int(start.value) - n)
    end.value = str(int(end.value) + n)
    # changing the values will trigger update_region_plot()


def update_gene(attr, old, new):

    gene_name = gene_search.value
    gene_match = genes[genes.name == gene_name]
    if len(gene_match != 0):
        gene_start = gene_match.start.iloc[0]
        gene_end = gene_match.end.iloc[0]
        # changing start and end will automatically trigger update_region_plot
        start.value = str(gene_start - 1500)
        end.value = str(gene_end + 1500)

    else:
        gene_search.value = 'No gene match'


def export_svg():
    p.output_backend ='svg'
    export_svgs(p, filename='bokeh_plot.svg')
    
############################ Read in data ######################################

# endo TSS expression
endo_tss_lb = pd.read_csv('../processed_data/endo_tss/lb/rLP5_Endo2_lb_expression_formatted_std.txt',
                           sep='\t')

# LB genomic shearing fragment pileup
frag_lb_plus = pd.read_csv('../processed_data/frag/lb/plus_frag_pileup.wig',
                            sep='\t', skiprows=1, names=['position', 'expression'])

frag_lb_minus = pd.read_csv('../processed_data/frag/lb/minus_frag_pileup.wig',
                            sep='\t', skiprows=1, names=['position', 'expression'])

# M9 minimal genomic shearing fragment pileup
frag_m9_plus = pd.read_csv('../processed_data/frag/m9/plus_frag_pileup_M9.wig',
                            sep='\t', skiprows=1, names=['position', 'expression'])

frag_m9_minus = pd.read_csv('../processed_data/frag/m9/minus_frag_pileup_M9.wig',
                            sep='\t', skiprows=1, names=['position', 'expression'])


# format frag pileup data
frag_plus_pileup = frag_lb_plus.merge(frag_m9_plus, on='position', how='outer', suffixes=['_lb', '_m9'])
frag_plus_pileup.columns = ['position', 'LB', 'M9']
frag_minus_pileup = frag_lb_minus.merge(frag_m9_minus, on='position', how='outer', suffixes=['_lb', '_m9'])
frag_minus_pileup.columns = ['position', 'LB', 'M9']

frag_plus_pileup = pd.melt(frag_plus_pileup, id_vars=['position'], value_vars=['LB', 'M9'],
             var_name='condition', value_name='expression')
frag_minus_pileup = pd.melt(frag_minus_pileup, id_vars=['position'], value_vars=['LB', 'M9'],
             var_name='condition', value_name='expression')

# normalize and scale
frag_plus_pileup_norm = deepcopy(frag_plus_pileup)
frag_plus_pileup_norm.expression = min_max_scale(frag_plus_pileup_norm.expression, NEW_MIN, NEW_MAX,
    median_norm = np.median(frag_plus_pileup_norm.expression))
frag_minus_pileup_norm = deepcopy(frag_minus_pileup)
frag_minus_pileup_norm.expression = min_max_scale(frag_minus_pileup_norm.expression, NEW_MIN, NEW_MAX,
    median_norm = np.median(frag_minus_pileup_norm.expression))


# read in M9 sequencing data
rna_plus = pd.read_csv('B6_M9_2.Forward.wig',
    sep='\t', skiprows=1, names=['position', 'raw_expression'])
rna_minus = pd.read_csv('B6_M9_2.Reverse.pos_values.wig',
    sep='\t', skiprows=1, names=['position', 'raw_expression'])
rna_plus.raw_expression.describe()

# normalize and scale
max_val = 1000
rna_plus = rna_plus.assign(expression = min_max_scale(rna_plus.raw_expression, NEW_MIN, NEW_MAX, max_value=max_val,
    median_norm=np.median(rna_plus.raw_expression)))
rna_minus = rna_minus.assign(expression = min_max_scale(rna_minus.raw_expression, NEW_MIN, NEW_MAX, max_value=max_val),
    median_norm=np.median(rna_minus.raw_expression))

# add condition
rna_plus['condition'] = ['M9'] * len(rna_plus)
rna_minus['condition'] = ['M9'] * len(rna_minus)


# read in gene annotation
genes = pd.read_csv('U00096.2_genes_clean.bed', sep = '\t', header = None,
                      names=['chrom', 'start', 'end', 'name', 'score', 'strand', 
                            'thick_start', 'thick_end', 'item_rgb', 'block_count',
                            'block_sizes', 'block_start'])
# drop unnecessary columns and simplify name
genes = genes[['start', 'end', 'name', 'strand']]
genes.name = [x.split(':')[1] for x in genes.name.tolist()]

# create colors for each condition and strand
colors = pd.DataFrame({'condition' : ['LB', 'LB', 'M9', 'M9'],
                      'strand' : ['+', '-', '+', '-'],
                      'color' : ['#8daec9', '#edac80', '#528ecb', '#ef8137']})


############################ Define widgets ####################################

condition_selection = CheckboxGroup(labels=['LB', 'M9'], active = [1])
# link change in selected buttons to update function
condition_selection.on_change('active', update_region_plot)

# # RangeSlider to change positions of region
# position_selection = RangeSlider(start = 2238500 - 300, end = 2238500 + 300,
#                           value = (2238500 - 300, 2238500 + 300),
#                           step = 100, title = 'genomic position')

# TextInput to define start and end position
# lacZ
# start = TextInput(value='362455', title='genome U00096.2 start position')
# end = TextInput(value='365529', title='genome end position')
start = TextInput(value='1232000', title='genome (U00096.2) start position')
end = TextInput(value='1238000', title='genome end position')

# update plot when value is changed
# position_selection.on_change('value', update_region_plot)
start.on_change('value', update_region_plot)
end.on_change('value', update_region_plot)

# toggle gene annotation on and off
gene_button = Select(
    title='RegulonDB gene annotation (U00096.2)',
    options=['show', 'hide'],
    value='show')
gene_button.on_change('value', toggle_genes)

# toggle TSS
tss_button = Select(
    title='Transcription Start Site (TSS) condition',
    options=['show', 'hide'],
    value='show')
tss_button.on_change('value', toggle_tss)

# toggle RNA
rna_button = Select(
    title='RNA-seq (M9)',
    options=['show', 'hide'],
    value='show')
rna_button.on_change('value', toggle_rna)

# add 500bp to both sides
add_button = Button(label='Add 500bp to both ends', button_type='primary')
add_button.on_click(add_bp)

# update axis range on click
axis_button = Button(label='Update Axis', button_type='primary')
axis_button.on_click(update_axis)

# search by gene name
gene_search = TextInput(value='Enter your favorite gene!', title='Search by gene name (case sensitive)')
gene_search.on_change('value', update_gene)

# export to SVG
export_button = Button(label='Export plot to SVG', button_type='primary')
export_button.on_click(export_svg)

############################## Initialize ######################################

# find initial conditions and position
initial_conditions = [condition_selection.labels[i] for i in condition_selection.active]
# split comma separated string, remove any leading or trialing whitespaces, convert to int
# initial_start, initial_end = map(int, map(unicode.strip, position_selection.value.split(',')))
initial_start = int(start.value)
initial_end = int(end.value)

src = make_pileup_dataset(initial_conditions,
    initial_start,
    initial_end,
    colors)

src_gene = make_region_genes(genes, initial_start, initial_end)
src_tss = make_tss_arrow(endo_tss_lb, initial_start, initial_end)
src_rna = make_rna_line(initial_start, initial_end, initial_conditions)

p = make_region_plot(src, src_gene, src_tss, src_rna)



################################# Layout ######################################

# combine all elements onto one page by defining layout

# add "figure" text
general_text = Div(text='''
    <b>General:</b> The data shown here is detailed in the paper "Genome-wide function characterization of
    E. coli promoters". In brief, we developed a massively parallel reporter assay that allows us to clone a
    library of DNA sequences upstream of a GFP reporter, which in turn drives 
    expression of a 20bp barcode. We can use next-generation sequencing to quantitatively measure activity of the
    entire library at the same time.''')

norm_text = Div(text='''
    <b>Normalization:</b> There are several normalization steps for visualization purposes to display 
    two different datasets (promoter activity and RNA-seq) on the same axis on the left hand side. 
    First, we remove values below the median by first subtracting the median from all values, and setting negative values to zero. 
    Next, we normalize and scale the data between [1, 1000]. Finally, we log-transform the data to 
    highlight small differences in activity. ''')

frag_text = Div(text='''
    <b>Promoter activity:</b> We created a library by randomly shearing the genome into small fragments (50-500bp) and cloned it
    into our MPRA, measuring the promoter activity of each fragment. We tested this 
    library in two different conditions (LB and M9). We condensed the data to show average activity
    at each nucleotide.''')

rna_text = Div(text='''
    <b>RNA-sequencing:</b> We performed genome-wide sequencing of E. coli in M9.''')

tss_text = Div(text='''
    <b> TSS library:</b> We synthesized a library of ~18K reported TSSs and measured their promoter activity.
    We synthesized 150bp surrounding the TSS (-120 to +30) in LB. We classified TSSs into active
    and inactive based on a set of 500 negative controls (150bp genomic sequence > 200bp from TSS). We set the active
    threshold at 3 standard deviations greater than the median negative control.
    ''')


# put controls in single element
# controls = WidgetBox(condition_selection, position_selection)
# controls_row1 = WidgetBox(general_text, norm_text, frag_text, condition_selection, 
    # rna_text, rna_button, axis_button)
# controls_row2 = WidgetBox(tss_text, tss_button, gene_button, start, end, axis_button)

# create row layout
# layout = gridplot([controls_row1, p, WidgetBox(tss_text, tss_button), WidgetBox(gene_button, start, end)], ncols=2)

controls_group1 = WidgetBox(general_text, norm_text, frag_text, condition_selection, rna_text, rna_button, axis_button)
controls_group2 = WidgetBox(gene_button, start, end, add_button, gene_search, export_button)
layout = row(controls_group1, p, controls_group2)

# make tab with layout
tab = Panel(child=layout, title = 'Genomic fragment shearing pileup')
tabs = Tabs(tabs=[tab])

# add to current document (displays plot)
curdoc().add_root(tabs)

