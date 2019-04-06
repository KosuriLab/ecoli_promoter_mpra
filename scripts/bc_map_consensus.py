"""
Simple barcode mapping utility.
Nathan Lubock

Outputs
"""

# Ensure Python 2/3 compatibility
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import itertools
import argparse
import multiprocessing
import sys
import csv
import os
import shlex
import subprocess
import time
import re
from signal import signal, SIGPIPE, SIG_DFL
from collections import Counter, defaultdict

# External Deps
import Levenshtein

# catch broken pipe errors to allow ex) python bc-map.py ... | head
# see: https://stackoverflow.com/a/30091579
signal(SIGPIPE, SIG_DFL)

#===============================================================================
# recepies/helper functions

def all_equal(iterable):
    "Returns True if all the elements are equal to each other"
    g = itertools.groupby(iterable)
    return next(g, True) and not next(g, False)

#===============================================================================
# I/O Functions

def fastq_reader(fastq):
    """
    Read in a fastq file laizy and return a generator of the sequence

    Parameters:
    -----------
    fastq :: FileType
        opened file to read

    Yields:
    -------
    generator :: str
        read portion of the fastq file

    Requires:
    ---------
    itertools
    """
    fourth = itertools.islice(fastq, 1, None, 4)
    for seq in fourth:
        yield seq.strip()

#-------------------------------------------------------------------------------

def fasta_reader(fasta):
    """
    Read in a fasta file lazily and return a generator of the name and sequence

    Parameters:
    -----------
    fasta :: FileType
        opened file

    Yields:
    -------
    generator :: (name, seq)
        name :: str
            Name of the read taken from the fasta file
        read :: str
            Sequence taken from the fasta file

    Requires:
    ---------
    itertools

    Example:
    --------
    itertools.groupby takes a key function and groups all items into a list
    until that key changes. We can key on lines beginning with >, then grab
    every line until the next record in the fasta. This makes our method robust
    to some fasta formats that have forced line breaks at given characters.

    foo = '>ABC>DEF>GHI'
    [(k, list(g)) for k,g in itertools.groupby(foo, lambda x: x == '>')]
    --> [(True, ['>']), (False, ['A', 'B', 'C']), (True, ['>']), ... ]

    Note:
    -----
    Adapted from: https://www.biostars.org/p/710/#1412
    """
    # ditch the boolean (x[0]) and just keep the header/seq grouping
    fa_iter = (x[1] for x in itertools.groupby(fasta, lambda line: line[0] == ">"))
    for header in fa_iter:
        # drop the ">"
        name = next(header)[1:].strip()
        # join all sequence lines to one by iterating until the next group.
        read = "".join(s.strip() for s in next(fa_iter))
        yield name, read

#===============================================================================
# Processing Functions

def get_idx(seq, bc_start, bc_len, var_start=None, var_len=None):
    """
    Get the slices corresponding to barcode and variant

    Parameters:
    -----------
    seq :: str
        Sequence to slice
    bc_start :: int
        Starting location of barcode (can be negatively indexed)
    bc_len :: int

    Optional arguments below are useful if there is junk sequence between
    variant and barcode
    var_start :: int
        Starting location of variant
    var_len :: int

    Returns:
    --------
    bc_slice :: slice
        location of the barcode
    var_slice :: slice
        location of the variant
    """
    if bc_start <= 0:
        bc_slice = slice(len(seq) + bc_start, len(seq) + bc_start + bc_len)
        if var_start is not None:
            # compensate for 0-based indexing
            var_slice = slice(var_start - 1, var_len)
        else:
            var_slice = slice(bc_start)
    else:
        # compensate for 0-based indexing
        bc_slice = slice(bc_start - 1, bc_start + bc_len - 1)
        var_slice = slice(bc_start + bc_len - 1, len(seq))
    return (bc_slice, var_slice)

#-------------------------------------------------------------------------------

def merge_reads(counter):
    """
    Merge list of strings into a consensus sequnce. Pad Ns for longer sequences.

    Parameters:
    -----------
    counter :: Counter
        Counter of strings corresponding to variants

    Returns:
    --------
    out_seq :: str
        Merged read

    Requires:
    ---------
    itertools
    all_equal -> from the itertools recipies

    Note:
    -----
    Assumes all reads are aligned at the left-most base (e.g. no indels at the
    very beginning of the sequence. Briefly, we transpose the sequences, grab
    the most common base at each position, and then trim the resulting
    consensus to the most common length (for >2 reads). For two or more reads,
    we will call an N at any bases without a consensus.
    """
    # one read
    if len(counter) == 1:
        return counter.keys()[0]

    # elements decompresses counter into iterable
    reads = counter.elements()
    trans = itertools.izip_longest(*reads, fillvalue='N')

    # eactly two reads (reads stored as counter vals)
    if sum(counter.itervalues()) == 2:
        consensus = (x[0] if x[0] == x[1] else 'N' for x in trans)
        return ''.join(consensus)
    else:
        trans_counter = (Counter(x) for x in trans)
        # return an N if the there are equal number of more than one nucleotide
        consensus = ('N' if all_equal(x.values()) and len(x.values()) > 1
                else x.most_common()[0][0] for x in trans_counter)
        raw_seq = ''.join(consensus)

        # get the lens to trim
        lens = itertools.chain.from_iterable(
                itertools.repeat(len(k), v) for k, v in counter.iteritems())
        common_len = Counter(lens).most_common()[0][0]
        return raw_seq[:common_len]

#-------------------------------------------------------------------------------

def mismatch(word, num_mismatches):
    """
    Generate all mismatches of a sequence upto a specified distance

    Parameters:
    -----------
    word :: str
    num_mismatches :: int

    Yields:
    -------
    _ :: str
        generator of sequences at num_mismatches from word

    Requires:
    ---------
    itertools

    Note:
    -----
    Boosted from: http://stackoverflow.com/a/11679968
    Briefly, generate all positions that will be edited. Insert a list of
    substitutions at those positions, then generate all possible combinations
    of that list.
    """
    letters = 'ATGC'
    for dist in range(1, num_mismatches + 1):
        for locs in itertools.combinations(range(len(word)), dist):
            this_word = [[char] for char in word]
            for loc in locs:
                orig_char = word[loc]
                this_word[loc] = [l for l in letters if l != orig_char]
            for poss in itertools.product(*this_word):
                yield ''.join(poss)

#-------------------------------------------------------------------------------

def bc_overlap(bcmap, dist):
    """
    Check if there are any sequences within a given distance to the input

    Parameters:
    -----------
    bcmap :: dict
        dict of sequences to check
    dist :: int
        distance to generate variants

    Returns:
    --------
    collector :: set
        Set of sequences from bcmap that are within dist of eachother

    Requires:
    ---------
    mismatch

    Note:
    -----
    Assumes that checking len(mismatch(seq, dist)) x len(bcmap) is less
    opperations than checking len(bcmap) x len(bcmap). Also uses a dict for
    constant time lookups.
    Briefly, iterate through every seq in bcmap, generate its variants, check
    to see if any of those variants are in bcmap, and then add them to a
    collector.
    """
    collector = set()
    for bc in bcmap:
        bc_var = mismatch(bc, dist)
        for var in bc_var:
            if var in bcmap and var not in collector:
                # assumes implicitly that bc -> var and var -> bc will be added
                collector.add(var)
    return collector


#-------------------------------------------------------------------------------

def intra_bc_contam(counter, nreads, cutoff):
    """
    Check if there is any sequence within a barcode that has >== nreads that are
    further than the most common read by some cuttoff

    Parameters:
    -----------
    counter :: Counter
        counter of reads and counts
    nreads :: int
        number of reads for ancillary sequence
    cutoff :: int
        Levenshtein distance cutoff
    Returns:
    --------
    bool
        True - sequence with >= nreads >= cutoff from most common seq
    Requires:
    ---------
    collections.Counter
    Levenshtein - https://pypi.python.org/pypi/python-Levenshtein/0.12.0
    """
    uniqd = counter.most_common()
    ancillary = [x[0] for x in uniqd[1:] if x[1] >= nreads]
    if len(uniqd) > 1 and ancillary:
        dists = (Levenshtein.distance(uniqd[0][0], x) for x in ancillary)
        if all(d < cutoff for d in dists):
            return False
        else:
            return True
    else:
        return False

#-------------------------------------------------------------------------------

def greater_diff(counter, cutoff):
    """
    Compare the most common item of a counter to everything else in it. If the
    difference is greater than some cutoff, return true

    Parameters:
    -----------
     counter :: Counter
        collection of numerics
    cutoff :: int
        difference cutoff (return false if abs(diff) >= cutoff)

    Returns:
    --------
    bool

    Requires:
    ---------
    collections.Counter
    """
    uniqd = counter.most_common()
    if len(uniqd) > 1:
        dists = [x[0] for x in uniqd]
        if any(abs(dists[0] - x) >= cutoff for x in dists):
            return True
        else:
            return False
    else:
        return False

#===============================================================================

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Collapse reads into a barcode-variant map. Barcodes are \
                filtered in two steps. First, we require a minimum number of \
                reads for every barcode. Second, we use BBMap to align the \
                reads associated with a given barcode. If the left-most \
                aligning base is  >= -d away from the most common read, then \
                we drop that barcode. Note that we expect the leftmost \
                position of all reads to bin into distinct values based on our \
                specific library design.')
    parser.add_argument('infile',
                        type=argparse.FileType('r'),
                        default=sys.stdin,
                        nargs='?',
                        help='path to a *.fastq file of the reads (or stdin if \
                                none)')
    parser.add_argument('ref',
                        type=str,
                        help='path to a fasta file of the reference for the variants')
    parser.add_argument('--rna_bcs',
                        type=argparse.FileType('r'),
                        help='path to a txt file of the rna-seq barcodes to keep.\
                        Optional, not necessary for initial mapping.')
    parser.add_argument('-v',
                        '--verbose',
                        dest='verbose',
                        action='store_true',
                        help='Verbose logging')
    parser.add_argument('-j',
                        '--proc',
                        dest='proc',
                        type=int,
                        default=1,
                        metavar='N',
                        choices=range(1, multiprocessing.cpu_count()+1),
                        help='number of processors (default=1, max={})'.format(multiprocessing.cpu_count()))
    parser.add_argument('-s',
                        '--bc-start',
                        dest='bc_start',
                        type=int,
                        default=1,
                        metavar='N',
                        help='starting position of barcode (1-indexed). Can \
                                also be negative for relative indexing from \
                                the end of a read (default=1)')
    parser.add_argument('-l',
                        '--bc-length',
                        dest='bc_length',
                        type=int,
                        default=20,
                        metavar='N',
                        help='length of barcode (default=20)')
    parser.add_argument('--var_start',
                        dest='var_start',
                        type=int,
                        help='starting position of variant (1-indexed). Provide if space between variant and barcode')
    parser.add_argument('--var_length',
                        dest='var_length',
                        type=int,
                        help='length of variant. Provide if space between variant and barcode')
    parser.add_argument('-m',
                        '--min-reads',
                        dest='min_reads',
                        type=int,
                        default=5,
                        metavar='N',
                        help='min consensus reads for each barcode (default=5)')
    parser.add_argument('-d',
                        '--start-dist',
                        dest='start_dist',
                        type=int,
                        metavar='N',
                        default=5,
                        choices=range(1,100),
                        help='ensure start of read for each barcode is within \
                                this distance (default 5)')
    parser.add_argument('-x',
                        '--bbmap-procs',
                        dest='bb_proc',
                        type=int,
                        default=1,
                        metavar='N',
                        choices=range(1, multiprocessing.cpu_count()+1),
                        help='number of processors for starcode portion \
                                (default=1, max={})'.format(multiprocessing.cpu_count()))
    parser.add_argument('-t',
                        '--trunc-len',
                        dest='trunc_len',
                        type=int,
                        metavar='N',
                        default=30,
                        choices=range(1,1000),
                        help='Drop barcode if any read length is +- N from the \
                                most common read (default=30).')
    parser.add_argument('-r',
                        '--contam-reads',
                        dest='contam_reads',
                        type=int,
                        metavar='N',
                        default=2,
                        choices=range(1,100),
                        help='Minimum number of reads a potential contaminated \
                                sequence must have for intra-barcode filtering \
                                (default=2).')
    parser.add_argument('-c',
                        '--contam-dist',
                        dest='contam_dist',
                        type=int,
                        metavar='N',
                        default=4,
                        choices=range(1,100),
                        help='distance cutoff for intra-barcode filtering \
                                (default=4)')
    parser.add_argument('-b',
                        '--bad-bcs',
                        dest='bad_bcs',
                        type=str,
                        metavar='foo',
                        help='output bad barcodes to this file')
    args = parser.parse_args()

    #---------------------------------------------------------------------------
    # DATA LOADING
    #---------------------------------------------------------------------------

    if args.verbose:
        timer = []
        timer.append(time.clock())
        print('Mapping BCs...', file=sys.stderr)

    # load the RNA-seq barcodes
    if args.rna_bcs:
        rna_bcs = set(line.rstrip() for line in args.rna_bcs)
    # if not rna_bcs:
    #     raise ValueError('RNA-seq barcode list must not be empty!')

    if args.verbose and args.rna_bcs:
            timer.append(time.clock())
            print('Loaded {} RNA-seq barcodes in {:.2f} seconds'.format(
                len(rna_bcs), timer[-1] - timer[-2]),
                file=sys.stderr)

    # only take barcodes that are present in the RNA-seq data
    bcmap = defaultdict(list)
    for seq in fastq_reader(args.infile):
        if args.var_start is not None:
            bc_slice, var_slice = get_idx(seq, args.bc_start, args.bc_length,
                args.var_start, args.var_length)
        else:
            bc_slice, var_slice = get_idx(seq, args.bc_start, args.bc_length)
        barcode = seq[bc_slice]
        if args.rna_bcs:
            if barcode in rna_bcs:
                bcmap[barcode].append(seq[var_slice])
        else:
            bcmap[barcode].append(seq[var_slice])

    # Compress map into counter after loading
    bcmap = {k:Counter(v) for k,v in bcmap.iteritems()}

    if args.verbose:
        timer.append(time.clock())
        print('Mapped {} barcodes in {:.2f} seconds'.format(
            len(bcmap.keys()), timer[-1] - timer[-2]),
            file=sys.stderr)
        print('Filtering BCs without a consensus of < {} reads'.format(args.min_reads),
              file=sys.stderr)

    #--------------------------------------------------------------------------
    # OUTPUT INTERMEDIATE MAP
    #--------------------------------------------------------------------------

    # bcmap = {k:Counter(v) for k,v in bcmap.iteritems() if sum(v.values()) > 2}
    # with open('big-dist-bc.txt', 'w') as f:
    #     print('bc\treads\tdist\tseq', file=f)
    #     for barcode,counter in bcmap.iteritems():
    #         uniqd = counter.most_common()
    #         if len(uniqd) > 1:
    #             dists = [0]
    #             seqs, reads = zip(*uniqd)
    #             dists.extend(Levenshtein.distance(seqs[0], x) for x in seqs[1:])
    #             if True in [d >= 50 for d in dists]:
    #                 for r,d,s in itertools.izip(reads, dists, seqs):
    #                     print('{}\t{}\t{}\t{}'.format(barcode,r,d,s), file=f)
    #
    # with open('bc-map.no-join.txt', 'w') as f:
    #     for k,v in bcmap.iteritems():
    #         for seq, read in v.iteritems():
    #             print('{}\t{}\t{}'.format(k, read, seq), file = f)

    #---------------------------------------------------------------------------
    # READ FILTER - force barcodes to have min amount of reads
    #---------------------------------------------------------------------------

    count = 0
    bad_bcs = []
    for k,v in bcmap.iteritems():
        # most_common(1) -> [(seq, reads)] for most common item in Counter
        # if v.most_common(1)[0][1] < args.min_reads:
        if sum(v.itervalues()) < args.min_reads:
            count += 1
            bad_bcs.extend((k, reads, 'min_reads', seq) for seq, reads in v.iteritems())

    # bcmap = {k:v for k,v in bcmap.iteritems() if v.most_common(1)[0][1] >= args.min_reads}
    bcmap = {k:v for k,v in bcmap.iteritems() if sum(v.itervalues()) >= args.min_reads}

    if args.verbose:
        timer.append(time.clock())
        print('Found {} BCs with < {} reads in {:.2f} secs'.format(
            count, args.min_reads, timer[-1] - timer[-2]),
            file=sys.stderr)
        print('Using BBMap to find contaminated barcodes at {} away'.format(
            args.start_dist),
            file=sys.stderr)


    #---------------------------------------------------------------------------
    # BBMAP - filter barcodes that map to different positions in ADRB2
    #---------------------------------------------------------------------------

    DEVNULL = open(os.devnull, 'w')
    command = shlex.split('bbmap.sh in=stdin.fasta fastareadlen=1000 out=stdout.sam \
            maxindel=200 nodisk=t noheader=t ref={} t={}'.format(
                args.ref, args.bb_proc))
    bbmap = subprocess.Popen(command,
                         stdout=subprocess.PIPE,
                         stdin=subprocess.PIPE,
                         stderr=DEVNULL)

    # subprocess.communicate() loads everything into memory instead of
    # streaming, but it avoids blocking when unix pipe fills up
    out_str = []
    for bc, counter in bcmap.iteritems():
        for seq, reads in counter.iteritems():
            out_str.append('>{}_{}\n{}'.format(bc, reads, seq))
    sam_blob = bbmap.communicate('\n'.join(out_str))[0]

    # use re.finditer to lazily split tabs since string is so big
    new_line = re.compile("[^\n]+")
    def lazy_split(string, pattern):
        return (x.group(0) for x in re.finditer(pattern, string))

    align_dist = defaultdict(list)
    for line in lazy_split(sam_blob, new_line):
        split = line.split('\t')
        bc, reads = split[0].split('_')
        align_dist[bc].append(itertools.repeat(int(split[3]), int(reads)))

    # link alignment position back to how many times that sequence showed up
    align_dist = {k:Counter(itertools.chain.from_iterable(v)) for k,v in align_dist.iteritems()}

    if args.verbose:
        timer.append(time.clock())
        print('Mapped barcodes in {:.2f} secs'.format(timer[-1] - timer[-2]),
            file=sys.stderr)

    #---------------------------------------------------------------------------
    # BARCODE CONTAMINATION - check for 3 cases:
    #   1) reads that align to different chunks
    #   2) truncations
    #   3) variants within the same chunk
    #---------------------------------------------------------------------------

    # process BBMap output for barcodes that align too far away
    contam_bcs = set(k for k,v in align_dist.iteritems() if greater_diff(v, args.start_dist))
    for bc in contam_bcs:
        for seq, reads in bcmap[bc].iteritems():
            bad_bcs.append((bc, reads, 'align', seq))
    bcmap = {k:v for k,v in bcmap.iteritems() if k not in contam_bcs}

    if args.verbose:
        timer.append(time.clock())
        print('Found {} barcodes that map too far away in {:.2f} secs'.format(
            len(contam_bcs), timer[-1] - timer[-2]),
            file=sys.stderr)

    #---------------------------------------------------------------------------
    # check for truncations by comparing to most common length
    contam_bcs = set(k for k,v in bcmap.iteritems() if
            greater_diff(Counter(len(x) for x in v.elements()), args.trunc_len))
    for bc in contam_bcs:
        for seq, reads in bcmap[bc].iteritems():
            bad_bcs.append((bc, reads, 'trunc', seq))
    bcmap = {k:v for k,v in bcmap.iteritems() if k not in contam_bcs}

    if args.verbose:
        timer.append(time.clock())
        print('Found {} truncated barcodes in {:.2f} secs'.format(
            len(contam_bcs), timer[-1] - timer[-2]),
            file=sys.stderr)

    #---------------------------------------------------------------------------
    # check for variants close to most common sequence in map
    contam_bcs = set(k for k,v in bcmap.iteritems() if
            intra_bc_contam(v, args.contam_reads, args.contam_dist))
    for bc in contam_bcs:
        for seq, reads in bcmap[bc].iteritems():
            bad_bcs.append((bc, reads, 'dist', seq))
    bcmap = {k:v for k,v in bcmap.iteritems() if k not in contam_bcs}

    if args.verbose:
        timer.append(time.clock())
        print('Found {} barcodes with at least one sequence at dist {} with at least {} reads in {:.2f} secs'.format(
            len(contam_bcs), args.contam_dist, args.contam_reads,
            timer[-1] - timer[-2]),
            file=sys.stderr)

    #---------------------------------------------------------------------------
    # READ MERGING
    #---------------------------------------------------------------------------

    def dummy(pack):
        """since we cant pickle lambdas define a dummy function to map over the tups"""
        key, value = pack
        num_reads = sum(value.itervalues())
        return (key, merge_reads(value), num_reads)

    # set up pool and distribute
    pool = multiprocessing.Pool(args.proc, maxtasksperchild=1000)
    writer = csv.writer(sys.stdout, lineterminator='\n')

    # ideally results are written as they come in from pool, unsure if true...
    consensus = pool.imap_unordered(dummy, bcmap.iteritems(), chunksize=10000)
    writer.writerows(consensus)

    pool.close()
    pool.join()

    if args.verbose:
        timer.append(time.clock())
        print('Merged reads in {:.2f} seconds'.format(timer[-1] - timer[-2]),
                file=sys.stderr)
        print('Total BCs remaining after filters {}'.format(len(bcmap)),
              file=sys.stderr)
        print('Total time {:.2f} mins'.format((time.clock() - timer[0])/60),
                file=sys.stderr)

    if args.bad_bcs is not None:
        write_start = time.clock()
        with open(args.bad_bcs, 'w') as f:
            writer = csv.writer(f, lineterminator='\n')
            writer.writerow(['bc', 'count', 'reason', 'var'])
            writer.writerows(bad_bcs)
        if args.verbose:
            print('Wrote {} barcodes in {:.2f} secs'.format(
                len(set(x[0] for x in bad_bcs)), time.clock() - write_start),
                file=sys.stderr)
