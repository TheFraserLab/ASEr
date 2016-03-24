#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Calculate gene/transcript-level counts from SNP counts.

============================================================================

        AUTHOR: Carlo Artieri
    MAINTAINER: Michael Dacre <mike.dacre@gmail.com>
  ORGANIZATION: Stanford University
       LICENSE: MIT License, property of Stanford, use as you wish
       CREATED: 2015-03-19
 Last modified: 2016-03-23 01:01

   DESCRIPTION: Calculates gene/transcript level counts from the output of
                CountSNPASE.py

============================================================================
"""
###########
# MODULES #
###########
import sys              # Access to simple command-line arguments
import argparse         # Access to long command-line parsing
from collections import defaultdict
from bisect import bisect
import pandas as pd

# Us
from ASEr import run    # File handling utilities

##########################
# COMMAND-LINE ARGUMENTS #
##########################

EPILOG = """\
NOTE:    SNPs that overlap multiple features on the same strand (or counting from unstranded
    libraries) will be counted in EVERY feature that they overlap. It is important to
    filter the annotation to count features of interest!

Detailed description of inputs/outputs follows:

-p/--phasedsnps
    A tab-delimited file (custom BED format) with positions and haplotype
    information of the masked SNPs of interest as follows:

    [CHR]\t[0 POSITION]\t[1 POSITION]\t[NAME]\t[REF|ALT]

    The fourth column MUST contain the phased SNPs alleles.

    This file can be created using the create_phased_bed script.

-g/--gff
    The script accepts both GTF and GFF annotation files. This should be combined with
    the -i/--identifier option specifying the identifier in the info column (column 9)
    that will be used for grouping counts. For example, in a GTF 'gene_id' will group
    counts by gene with 'transcript_id' with group counts by transcript. In addition,
    the -t/--type option sets the feature type (column 3) from which to pull features
    typically you'd want to count from 'exon', but many annotations may use non-standard
    terms.

-m/--min
    This sets the minimum # of reads required to include a SNP in the calculation of the
    fraction of SNPs agreeing in allelic direction.

-w/--writephasedsnps
    If this is specified, then the program will output an additional output file named
    [OUTFILE].snp.txt with phased SNP-level ASE calls. This can be useful for checking
    SNP consistency across samples. See below for a description of the output.

-s/--stranded
    If the data come from a stranded library prep, then this option will only count reads
    mapped to the corresponding strand.

OUTPUT:

The output of the script is a tab-delimited text file set by -o/--outfile, which contains
the following columns:

FEATURE         Name of the counted feature
CHROMOSOME         Chromosome where feature is found
ORIENTATION         Orientation of feature (+/-)
START-STOP         Ultimate 5' and 3' 1-based start and stop positions
REFERENCE_COUNTS     Total reference allele counts across SNPS (or first allele in the REF|ALT phasing)
ALT_COUNTS         Total alternate allele counts across SNPs (or second allele in the REF|ALT phasing)
TOTAL_SNPS         The total number of SNPs overlapped by the feature
REF_BIASED         Number of REF biased SNPs passing the -m/--min threshold
ALT_BIASED         Number of ALT biased SNPs passing the -m/--min threshold
REF-ALT_RATIO         The proportion of SNPs agreeing in direction (0.5 - 1)
SNPS             A list of all SNPs overlapped by the feature separated by ';' and of the format:

    [1-based position],[REF_ALLELE]|[ALT_ALLELE],[REF_COUNTS]|[ALT_COUNTS];

If the -w/--writephasedsnps option has been set, it will produce a tab-delimited table with the
following columns:

CHROMOSOME         Chromosome where SNP is found
POSITION         1-based position
FEATURE         Feature in which SNP is found
ORIENTATION         Orientation of feature (if stranded only reads on this strand are counted)
REFERENCE_ALLELE     Reference base
ALTERNATE_ALLELE     Alternate base
REF_COUNTS         Reference base counts
ALT_COUNTS         Alternate base counts


"""


def read_snp_count_file(snp_file):
    """Return a position data structure from a SNP counts file.

    The SNP counts file come from CountSNPASE

    :returns: A dictionary with the format:
                position => {pos_dict[BASE] => positive strand count
                             neg_dict[BASE] => negative strand count

    """
    snp_counts_dict = {}
    with run.open_zipped(snp_file) as count_file:
        for line in count_file:
            if 'SUM_POS_READS' in line:
                continue

            line   = line.rstrip('\n')
            line_t = line.split('\t')

            pos = line_t[0] + '|' + str(line_t[1])

            pos_counts = line_t[2].split('|')
            neg_counts = line_t[3].split('|')
            bases = ['A', 'C', 'G', 'T']

            pos_dict = {}
            neg_dict = {}

            for i in range(len(pos_counts)):
                pos_dict[bases[i]] = pos_counts[i]
                neg_dict[bases[i]] = neg_counts[i]

            snp_counts_dict[pos] = [pos_dict, neg_dict]
    return snp_counts_dict


def read_snp_phasing_file(snp_file):
    """Return a dictionary of phased SNPs by genome position.

    :snp_file: A bed format file with phased SNPs, can be produces with the
               create_phased_bed script.
    :returns:  A dictionary of::
                    chromosome => [(pos, REF|ALT), ...]
    """
    snp_phase_dict = defaultdict(list)
    with run.open_zipped(snp_file) as snp_file:
        for line in snp_file:
            line   = line.rstrip('\n')
            line_t = line.split('\t')

            chrom = line_t[0]
            pos = int(line_t[2])
            if len(snp_phase_dict[chrom]) and pos < snp_phase_dict[chrom][-1][0]:
                raise ValueError("SNP file {} should be sorted by chromosome position"
                        .format(snp_file.name))
            snp_phase_dict[chrom].append((pos, line_t[3]))
    return snp_phase_dict


def main(argv=None):
    """Run as a script."""
    if not argv:
        argv = sys.argv[1:]

    parser = argparse.ArgumentParser(
        description='Takes the output of CountSNPASE.py and generates gene level ASE counts.',
        add_help=False, epilog=EPILOG, formatter_class=run.CustomFormatter)

    req = parser.add_argument_group('Required arguments:')
    req.add_argument('-c', '--snpcounts', action="store", dest="snpcounts",
                     help='SNP-level ASE counts from CountSNPASE.py ',
                     required=True, metavar='')
    req.add_argument('-p', '--phasedsnps', action="store", dest="phasedsnps",
                     help='BED file of phased SNPs', required=True, metavar='')
    req.add_argument('-g', '--gff', action="store", dest="gff",
                     help='GFF/GTF formatted annotation file', required=True,
                     metavar='')
    req.add_argument('-o', '--outfile', action="store", dest="outfile",
                     help='Gene-level ASE counts output', required=True,
                     metavar='')

    opt = parser.add_argument_group('Optional arguments:')
    opt.add_argument('-w', '--writephasedsnps', action="store_true", dest="write",
                     help='Write a phased SNP-level ASE output file [OUTFILE].snps.txt')
    opt.add_argument('-i', '--identifier', action="store", dest="id",
                     help='ID attribute in information column',
                     default='gene_id', metavar='')
    opt.add_argument('-t', '--type', action="store", dest="type",
                     help='Annotation feature type', default='exon', metavar='')
    opt.add_argument('-m', '--min', action="store", dest="min", type=int,
                     help='Min reads to calculate proportion ref/alt biased',
                     default=10)
    opt.add_argument('-s', '--stranded', action="store_true", dest="stranded",
                     help='Data are stranded? [Default: False]')

    opt.add_argument('-h', '--help', action="help",
                     help="Show this help message and exit")

    args = parser.parse_args()

    ##########
    # SCRIPT #
    ##########

    # Read in the SNP-level ASE counts
    snp_counts_dict = read_snp_count_file(args.snpcounts)

    # Read in the SNP phasing information
    snp_phase_dict = read_snp_phasing_file(args.phasedsnps)

    # Initialize variables for gene-level counts
    chromosome       = {}
    position         = {}
    ori              = {}
    total_ref        = {}
    total_alt        = {}
    total_snps       = {}
    ref_biased       = {}
    alt_biased       = {}
    snp_array        = {}
    phased_snp_array = []

    features = {}    # We have to store a list of all of the features somewhere

    # Now go through the GFF/GTF and concatenate the approporiate info
    gff_file = open(args.gff, 'r')
    with run.open_zipped(args.gff) as gff_file:
        for line in gff_file:
            line   = line.rstrip('\n')
            line_t = line.split('\t')

            if line_t[2] != args.type:
                continue

            line_t8 = line_t[8].split(';')

            # Check for custom identifier and determine filetype:
            if args.id + '=' in line_t[8]:    # GFF
                for i in line_t8:
                    i.strip()
                    if args.id + '=' in i:
                        if args.id in i:
                            i2 = i.split('=')
                            name = i2[1]

            elif args.id + ' ' in line_t[8]:  # GTF
                for i in line_t8:
                    i.strip()
                    if args.id + ' ' in line_t[8]:
                        if args.id in i:
                            i2 = i.split('"')
                            name = i2[1]

            else:
                sys.stderr.write(('ID attribute "{}" doesn\'t exist or GFF/GTF file ' +
                                  'not properly formatted.\n').format(args.id))
                sys.stderr.write('GFF info column format is: attribute=value;\n')
                sys.stderr.write('GTF info column format is: attribute "value";\n')
                sys.exit(1)

            features[name] = 1
            chrom = line_t[0]
            chromosome[name] = chrom

            ori[name] = line_t[6]
            orientation = line_t[6]

            if name not in position:
                position[name] = []
                position[name].append(int(line_t[3]))
                position[name].append(int(line_t[4]))
            else:
                position[name].append(int(line_t[3]))
                position[name].append(int(line_t[4]))

            # Now go through the positions overlapped by the annotation and add in SNPs if appropriate
            if chrom not in snp_phase_dict:
                continue
            lower = bisect(snp_phase_dict[chrom], (int(line_t[3]), ))
            upper = bisect(snp_phase_dict[chrom], (int(line_t[4])+1, ))
            for i, snp in snp_phase_dict[chrom][lower:upper]:
                pos = chrom + '|' + str(i)

                if pos in snp_counts_dict:

                    # Count SNPs
                    if name in total_snps:
                        total_snps[name] += 1
                    else:
                        total_snps[name] = 0
                        total_snps[name] += 1

                    # Get REF|ALT counts
                    # Parse the REF|ALT dict
                    refalt = snp.split('|')

                    # Get pos/neg counts
                    ref_pos_counts = snp_counts_dict[pos][0][refalt[0]]
                    ref_neg_counts = snp_counts_dict[pos][1][refalt[0]]
                    alt_pos_counts = snp_counts_dict[pos][0][refalt[1]]
                    alt_neg_counts = snp_counts_dict[pos][1][refalt[1]]

                    # If stranded add only appropriate strand counts
                    if args.stranded is True:
                        if orientation == '+':
                            if name in total_ref:
                                total_ref[name] += int(ref_pos_counts)
                            else:
                                total_ref[name] = 0
                                total_ref[name] += int(ref_pos_counts)

                            if name in total_alt:
                                total_alt[name] += int(alt_pos_counts)
                            else:
                                total_alt[name] = 0
                                total_alt[name] += int(alt_pos_counts)

                            if name not in ref_biased:
                                ref_biased[name] = 0
                            if name not in alt_biased:
                                alt_biased[name] = 0

                            # Determine if ref or alt biased
                            if ref_pos_counts + alt_pos_counts >= args.min:
                                if ref_pos_counts > alt_pos_counts:
                                    ref_biased[name] += 1

                                elif ref_pos_counts < alt_pos_counts:
                                    alt_biased[name] += 1

                            # Add it to the total SNP arrays
                            if name in snp_array:
                                snp_array[name].append(
                                    str(i) + ',' + str(snp) +
                                    ',' + str(ref_pos_counts) + '|' +
                                    str(alt_pos_counts))
                                phased_snp_array.append(
                                    str(chromosome[name]) + '\t' + str(i) +
                                    '\t' + name + '\t' + ori[name] + '\t' +
                                    str(refalt[0]) + '\t' + str(refalt[1]) +
                                    '\t' + str(ref_pos_counts) + '\t' +
                                    str(alt_pos_counts))
                            else:
                                snp_array[name] = []
                                snp_array[name].append(
                                    str(i) + ',' + str(snp) +
                                    ',' + str(ref_pos_counts) + '|' +
                                    str(alt_pos_counts))
                                phased_snp_array.append(
                                    str(chromosome[name]) + '\t' + str(i) +
                                    '\t' + name + '\t' + ori[name] + '\t' +
                                    str(refalt[0]) + '\t' + str(refalt[1]) +
                                    '\t' + str(ref_pos_counts) + '\t' +
                                    str(alt_pos_counts))

                        elif orientation == '-':
                            if name in total_ref:
                                total_ref[name] += int(ref_neg_counts)
                            else:
                                total_ref[name] = 0
                                total_ref[name] += int(ref_neg_counts)

                            if name in total_alt:
                                total_alt[name] += int(alt_neg_counts)
                            else:
                                total_alt[name] = 0
                                total_alt[name] += int(alt_neg_counts)

                            if name not in ref_biased:
                                ref_biased[name] = 0
                            if name not in alt_biased:
                                alt_biased[name] = 0

                            # Determine if ref or alt biased
                            if ref_neg_counts + alt_neg_counts >= args.min:
                                if ref_neg_counts > alt_neg_counts:
                                    ref_biased[name] += 1

                                elif ref_neg_counts < alt_neg_counts:
                                    alt_biased[name] += 1

                            # Add it to the total SNP array
                            if name in snp_array:
                                snp_array[name].append(
                                    str(i) + ',' + str(snp) +
                                    ',' + str(ref_neg_counts) + '|' +
                                    str(alt_neg_counts))
                                phased_snp_array.append(
                                    str(chromosome[name]) + '\t' + str(i) +
                                    '\t' + name + '\t' + ori[name] + '\t' +
                                    str(refalt[0]) + '\t' + str(refalt[1]) +
                                    '\t' + str(ref_neg_counts) + '\t' +
                                    str(alt_neg_counts))
                            else:
                                snp_array[name] = []
                                snp_array[name].append(
                                    str(i) + ',' + str(snp) +
                                    ',' + str(ref_neg_counts) + '|' +
                                    str(alt_neg_counts))
                                phased_snp_array.append(
                                    str(chromosome[name]) + '\t' + str(i) +
                                    '\t' + name + '\t' + ori[name] + '\t' +
                                    str(refalt[0]) + '\t' + str(refalt[1]) +
                                    '\t' + str(ref_neg_counts) + '\t' +
                                    str(alt_neg_counts))

                    else:
                        if name in total_ref:
                            total_ref[name] += int(ref_pos_counts)
                            total_ref[name] += int(ref_neg_counts)
                        else:
                            total_ref[name] = 0
                            total_ref[name] += int(ref_pos_counts)
                            total_ref[name] += int(ref_neg_counts)

                        if name in total_alt:
                            total_alt[name] += int(alt_pos_counts)
                            total_alt[name] += int(alt_neg_counts)
                        else:
                            total_alt[name] = 0
                            total_alt[name] += int(alt_pos_counts)
                            total_alt[name] += int(alt_neg_counts)

                        # Determine if ref or alt biased
                        if name not in ref_biased:
                            ref_biased[name] = 0

                        if name not in alt_biased:
                            alt_biased[name] = 0

                        tot_ref = int(ref_pos_counts) + int(ref_neg_counts)
                        tot_alt = int(alt_pos_counts) + int(alt_neg_counts)
                        if tot_ref + tot_alt >= args.min:
                            if tot_ref > tot_alt:
                                if name in ref_biased:
                                    ref_biased[name] += 1
                                else:
                                    ref_biased[name] = 0
                                    ref_biased[name] += 1

                                if name in alt_biased:
                                    pass
                                else:
                                    alt_biased[name] = 0

                            elif tot_ref < tot_alt:
                                if name in alt_biased:
                                    alt_biased[name] += 1
                                else:
                                    alt_biased[name] = 0
                                    alt_biased[name] += 1

                                if name in ref_biased:
                                    pass
                                else:
                                    ref_biased[name] = 0

                        # Add it to the total SNP array
                        if name in snp_array:
                            snp_array[name].append(
                                str(i) + ',' + str(snp) +
                                ',' + str(int(tot_ref)) + '|' +
                                str(int(tot_alt)))
                        else:
                            snp_array[name] = []
                            snp_array[name].append(
                                str(i) + ',' + str(snp) + ',' +
                                str(int(tot_ref)) + '|' + str(int(tot_alt)))

    # Print the output
    keys = sorted(list(features.keys()))

    with open(args.outfile, 'w') as outfile:
        posits = {key:'{}-{}'.format(min(position[key]), max(position[key])) for key in keys}
        snp_arrays = {key: ';'.join(snp_array[key]) for key in total_ref}
        out = pd.DataFrame({
            'CHROMOSOME': chromosome,
            'ORIENTATION': ori,
            'START-STOP': posits,
            'REFERENCE_COUNTS': total_ref,
            'ALT_COUNTS' : total_alt,
            'TOTAL_SNPS': total_snps,
            'REF_BIASED' : ref_biased,
            'ALT_BIASED' : alt_biased,
            'REF-ALT_RATIO': pd.np.nan,
            'SNPS' : snp_arrays},
            index=keys,
            columns=['CHROMOSOME', 'ORIENTATION', 'START-STOP',
                'REFERENCE_COUNTS', 'ALT_COUNTS', 'TOTAL_SNPS', 'REF_BIASED',
                'ALT_BIASED', 'REF-ALT_RATIO', 'SNPS'])
        out.index.name = 'FEATURE'
        out['REF-ALT_RATIO'] = (out.ix[:, ['REF_BIASED', 'ALT_BIASED']].max(axis=1)
                / out.ix[:, ['REF_BIASED', 'ALT_BIASED']].sum(axis=1))
        out.ix[(out.REF_BIASED == 0) & (out.ALT_BIASED == 0), 'REF-ALT_RATIO'] = 1
        out.to_csv(outfile, sep='\t', na_rep='NA', float_format='%.16g')




    if args.write is True:
        with open(args.outfile + '.snps.txt', 'w') as outfile:
            # Header
            outfile.write('CHROMOSOME\tPOSITION\tFEATURE\tORIENTATION\t' +
                          'REFERENCE_ALLELE\tALTERNATE_ALLELE\tREF_COUNTS\t' +
                          'ALT_COUNTS\n')

            for i in phased_snp_array:
                outfile.write(i + '\n')

if __name__ == '__main__' and '__file__' in globals():
    sys.exit(main())
