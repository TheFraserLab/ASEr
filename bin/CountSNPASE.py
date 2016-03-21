#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Count number of reads overlapping each SNP in a sam/bam file.

============================================================================

        AUTHOR: Carlo Artieri
    MAINTAINER: Michael D Dacre, mike.dacre@gmail.com
  ORGANIZATION: Stanford University
       LICENSE: MIT License, property of Stanford, use as you wish
       VERSION: 0.1
       CREATED: 2015-03-16
 Last modified: 2016-03-21 07:59

   DESCRIPTION: This script will take a BAM file mapped to a SNP-masked
                genome and count the number of reads overlapping each SNP.
                When reads overlap multiple SNPs, the read is randomly
                assigned to one of them.

============================================================================
"""
import os                 # Path manipulation
import sys                # Access to simple command-line arguments
import argparse           # Access to long command-line parsing
import datetime           # Access to calendar/clock functions
import re                 # Access to REGEX splitting
import math               # Access to math functions
import random             # Access to random number generation
import subprocess         # Access to direct command line in/out processing
import textwrap           # Add text block wrapping properties
from time import sleep    # Allow system pausing

# Us
from ASEr import logme    # Logging functions
from ASEr import run      # File handling functions
from ASEr import cluster  # Queue submission

# Logging
logme.MIN_LEVEL = 'info'  # Switch to 'debug' for more verbose, 'warn' for less

###############################################################################
#                          Command Line Description                           #
###############################################################################

epilog = """\
Detailed description of inputs/outputs follows:

-s/--snps
    A tab-delimited BED file with positions of masked SNPs of interest as follows:

    [CHR]    [0 POSITION]    [1 POSITION]

    Additional columns are ignored.

-r/--reads
    A SAM or BAM file containing all of the reads masked to the masked genome. The file
    shound have all duplicates removed and MUST be sorted by read name
    (i.e. samtools sort -n ).

-m/--mode
    The script can be run in two modes. In 'single' mode, the entire SNP counting is performed
    locally. In 'multi' mode, the read file will be split up by the number of specified jobs on
    the cluster. This is much faster for large SAM/BAM files.

OUTPUT:

The output of the script is a tab-delimited text file, [PREFIX]_SNP_COUNTS.txt, which contains the
following columns:

CHR\t\tChromosome where SNP is found
POSITION\t1-based position of SNP
POS_A|C|G|T\tCount of reads containing A|C|G|T bases at the SNP position on the POSITIVE strand
NEG_A|C|G|T\tCount of reads containing A|C|G|T bases at the SNP position on the NEGATIVE strand
SUM_POS_READS\tSum of all reads assigned to the SNP on POSITIVE strand
SUM_NEG_READS\tSum of all reads assigned to the SNP on NEGATIVE strand
SUM_READS\tSum of all reads assigned to the SNP

"""

#############
# FUNCTIONS #
#############


def fasta_to_dict(file):
    """Convert a FASTA file to a dictionary.

    keys = headers and values = sequence.
    """
    fasta_file = run.open_zipped(file, "r")    # Open the file for reading
    fasta_dict = {}

    for line in fasta_file:
        line = line.rstrip('\n')
        if re.match('^>', line):
            line_split = line.split(' ')
            header = line_split[0].translate(None, '>')
            fasta_dict[header] = ''
        else:
            fasta_dict[header] += line

    fasta_file.close()

    return fasta_dict


def split_CIGAR(cigar):
    """Split the CIGAR string and return two lists: types and values in order."""
    cig_types_tmp = re.split('[0-9]', cigar)
    cig_vals_tmp = re.split('[MIDNSHP\=X]', cigar)
    cig_types = []
    cig_vals = []

    for i in cig_types_tmp:
        if i != '':
            cig_types.append(i)

    for i in cig_vals_tmp:
        if i != '':
            cig_vals.append(i)

    return cig_types, cig_vals


def CIGAR_to_Genomic_Positions(cigar_types, cigar_vals, pos):
    """Use the CIGAR string to return a list of genomic coordinates.

    Coordinates correspond to the individual bases of the read to get SNP
    positions from the MD tag.
    """
    # Initialize the list of genomic positions
    genomic_positions = []
    curr_pos = pos

    for i in range(len(cigar_types)):
        # What are we going to do to each CIGAR str.
        if cigar_types[i] == 'N':
            curr_pos = int(curr_pos) + int(cigar_vals[i])
        elif cigar_types[i] == 'D':
            curr_pos = int(curr_pos) + int(cigar_vals[i])
        elif cigar_types[i] == 'M':
            genomic_positions = genomic_positions + list(range(int(curr_pos),
                                                               int(curr_pos)+int(cigar_vals[i])))
            curr_pos = int(curr_pos) + int(cigar_vals[i])
    return genomic_positions


class CustomFormatter(argparse.ArgumentDefaultsHelpFormatter,
                      argparse.RawDescriptionHelpFormatter):

    """Custom argparse formatting."""

    pass


def split_samfile(sam_file, splits, prefix='', path=''):
    """Take a sam file and split it splits number of times.

    :path:    Where to put the split files.
    :prefix:  A prefix for the outfile names.
    :returns: A tuple of job files.
    """
    # Determine how many reads will be in each split sam file.
    num_lines = os.popen('wc -l ' + sam_file + ' | awk \'{print $1}\'').read()
    num_reads = int(int(num_lines)/splits) + 1

    # Subset the SAM file into X number of jobs
    cnt      = 0
    currjob  = 1
    suffix   = '.split_sam_' + str(currjob).zfill(4)
    sam_path, sam_name = os.path.split(sam_file)
    run_file = os.path.join(path, prefix + sam_name + suffix)
    outfiles = [run_file]

    # Actually split the file
    with open(sam_file) as in_sam:
        sam_split = open(run_file, 'w')
        for line in in_sam:
            cnt += 1
            if cnt < num_reads:
                sam_split.write(line)
            elif cnt == num_reads:
                line_t = line.split('\t')

                # Check if next line is mate-pair. If so, don't split across files.
                line2 = next(in_sam)
                line2_t = line2.split('\t')

                if line_t[0] == line2_t[0]:
                    sam_split.write(line)
                    sam_split.write(line2)
                    sam_split.close()
                    currjob += 1
                    suffix = '.split_sam_' + str(currjob).zfill(4)
                    sam_split = open(prefix + sam_file + suffix, 'w')
                    cnt = 0
                else:
                    sam_split.write(line)
                    sam_split.close()
                    currjob += 1
                    suffix = '.split_sam_' + str(currjob).zfill(4)
                    run_file = os.path.join(path, prefix + sam_name + suffix)
                    sam_split = open(run_file, 'w')
                    outfiles.append(run_file)
                    sam_split.write(line2)
                    cnt = 0
        sam_split.close()
    return tuple(outfiles)


def main(argv=None):
    """Main script."""
    ##########################
    # COMMAND-LINE ARGUMENTS #
    ##########################

    # Get myself
    program_name = sys.argv[0]

    if not argv:
        argv = sys.argv[1:]

    parser  = argparse.ArgumentParser(
        description=__doc__,
        add_help=False, epilog=epilog, formatter_class=CustomFormatter)

    req = parser.add_argument_group('Required arguments')
    req.add_argument('-m', '--mode',
                     help='Operation mode', choices=['single', 'multi'],
                     required=True, metavar='mode')
    req.add_argument('-s', '--snps',
                     help='SNP BED file', required=True, metavar='<BED>')
    req.add_argument('-r', '--reads',
                     help='Mapped reads file [sam or bam]',
                     required=True, metavar='<[S/B]AM>')

    uni = parser.add_argument_group('Universal optional arguments')
    uni.add_argument('-p', '--prefix',
                     help='Prefix for temp files and output', default='TEST',
                     metavar='')
    uni.add_argument('-b', '--bam', action='store_true', dest='bam',
                     help='Mapped read file type is bam (auto-detected if *.bam)')
    uni.add_argument('-t', '--single', action='store_true', dest='single',
                     help='Mapped reads are single-end')
    uni.add_argument('-n', '--noclean', action='store_true',
                     help='Do not delete intermediate files (for debuging)')
    uni.add_argument('-h', '--help', action='help',
                     help='show this help message and exit')

    mult = parser.add_argument_group('Multi(plex) mode arguments')
    mult.add_argument('-j', '--jobs', type=int,
                      help='Divide into # of jobs', default=100, metavar='')
    mult.add_argument('-w', '--walltime',
                      help='Walltime for each job', default='3:00:00',
                      metavar='')
    mult.add_argument('-k', '--mem', dest='memory',
                      help='Memory for each job', default='5000MB', metavar='')
    mult.add_argument('--queue',
                      help='Queue to submit jobs to', default='batch',
                      metavar='')
    mult.add_argument('--cluster', help='Which cluster to use',
                      choices=['torque', 'slurm'], default='torque',
                      metavar='')

    single = parser.add_argument_group('Single mode arguments')
    single.add_argument('-f', '--suffix',
                        help='Suffix for multiplexing [set automatically]',
                        default='', metavar='')

    logging = parser.add_argument_group('Logging options')
    logging.add_argument('-q', '--quiet', action='store_true',
                         help="Quiet mode, only prints warnings.")
    logging.add_argument('-v', '--verbose', action='store_true',
                         help="Verbose mode, prints debug info too.")
    logging.add_argument('--logfile',
                         help='Logfile to write messages too, default is ' +
                         'STDERR')

    args = parser.parse_args()

    ###########################################################################
    #                            File Preparations                            #
    ###########################################################################

    # Take care of logging
    if args.logfile:
        logme.LOGFILE = args.logfile
    if args.quiet:
        logme.MIN_LEVEL = 'warn'
    elif args.verbose:
        logme.MIN_LEVEL = 'debug'

    # Initialize variables
    prefix = args.prefix + '_'
    wasbam = False

    # Make sure we can run ourselves
    if not run.is_exe(program_name):
        program_name = run.which(parser.prog)

    # Set the cluster type if we are in multi mode
    if args.mode == 'multi':
        cluster.QUEUE = args.cluster

    # Check if the read file is sam or bam
    file_check = args.reads.split('.')
    file_check[-1] = file_check[-1].lower()
    sam_path, sam_file = os.path.split(args.reads)

    if file_check[-1] == 'bam' or args.bam is True:
        sam_file = '.'.join(args.reads.split('.')[:-1]) + '.sam'
        wasbam = True
        if not os.path.isfile(sam_file):
            logme.log('Converting BAM to SAM ...')
            os.system('samtools view ' + args.reads + ' > ' + sam_file)
            logme.log('Done')
        else:
            logme.log('BAM to SAM conversion already complete, using SAM ' +
                      'file.')

    ##################
    # MULTIPLEX MODE #
    ##################

    # If we're running in multiplex mode
    if args.mode == 'multi':
        logme.log('Splitting sam file {} into {} files.'.format(sam_file,
                                                                args.jobs))
        reads_files = split_samfile(sam_file, args.jobs, prefix)
        logme.log('Splitting complete.')

        # Create PBS scripts and submit jobs to the cluster
        subnoclean = '--noclean' if args.noclean else ''
        job_ids = []
        logme.log('Submitting split files to cluster')
        for reads_file in reads_files:
            suffix = reads_file[-4:]

            command = ("python2 " + program_name + " --mode single --snps " +
                       args.snps + " --reads " + reads_file + " --suffix " +
                       suffix + " --prefix " + args.prefix)

            job_ids.append(cluster.submit(command, name=prefix + suffix,
                                          time=args.walltime, cores=1,
                                          mem=args.memory,
                                          partition=args.queue))
            sleep(2)    # Pause for two seconds to make sure job is properly submitted

        # Now wait and check for all jobs to complete every so long
        logme.log('Submission done, waiting for jobs to complete.')
        done = False
        while done is False:
            tot_done = 0
            for i in range(1, args.jobs+1):
                suffix = str(i).zfill(4)

                if os.path.isfile(prefix + suffix + '_done'):
                    tot_done += 1

            if tot_done == args.jobs:
                done = True

            sleep(10)

        logme.log('Jobs completed.')
        os.system('rm *_done')    # Remove 'done' files in case we want to run again.

        # Once the jobs are done, concatenate all of the counts into one file.
        # Initialize dictionaries

        tot_pos_counts = {}
        tot_neg_counts = {}
        tot_tot_counts = {}
        tot_sum_pos = {}
        tot_sum_neg = {}

        for i in range(1, args.jobs+1):
            suffix = str(i).zfill(4)
            in_counts = run.open_zipped(prefix + 'SNP_COUNTS_' + suffix, 'r')

            # Parse the line to add it to the total file
            for line in in_counts:
                line = line.rstrip('\n')
                line_t = line.split('\t')

                if 'CHR' in line:
                    continue

                pos = line_t[0] + '|' + line_t[1]

                pos_split = line_t[2].split('|')
                neg_split = line_t[3].split('|')

                if pos in tot_pos_counts or pos in tot_neg_counts or pos in tot_tot_counts:
                    for j in range(len(pos_split)):
                        tot_pos_counts[pos][j] += int(pos_split[j])
                        tot_neg_counts[pos][j] += int(neg_split[j])
                    tot_sum_pos[pos] += int(line_t[4])
                    tot_sum_neg[pos] += int(line_t[5])
                    tot_tot_counts[pos] += int(line_t[6])

                else:
                    tot_pos_counts[pos] = [0, 0, 0, 0]
                    tot_neg_counts[pos] = [0, 0, 0, 0]
                    tot_tot_counts[pos] = 0
                    tot_sum_pos[pos] = 0
                    tot_sum_neg[pos] = 0
                    for j in range(len(pos_split)):
                        tot_pos_counts[pos][j] += int(pos_split[j])
                        tot_neg_counts[pos][j] += int(neg_split[j])
                    tot_sum_pos[pos] += int(line_t[4])
                    tot_sum_neg[pos] += int(line_t[5])
                    tot_tot_counts[pos] += int(line_t[6])

            in_counts.close()

        # Write out the final concatenated file
        final_counts = run.open_zipped(prefix + 'SNP_COUNTS.txt', 'w')
        final_counts.write('CHR\tPOSITION\tPOS_A|C|G|T\tNEG_A|C|G|T\tSUM_POS_READS\tSUM_NEG_READS\tSUM_READS\n')

        keys = sorted(tot_pos_counts.keys())

        for key in keys:
            pos = key.split('|')
            pos_fix = [str(x) for x in tot_pos_counts[key]]
            neg_fix = [str(x) for x in tot_neg_counts[key]]
            pos_out = '|'.join(pos_fix)
            neg_out = '|'.join(neg_fix)
            final_counts.write(str(pos[0]) + '\t' + str(pos[1]) + '\t' + pos_out + '\t' + neg_out + '\t' + str(tot_sum_pos[key]) + '\t' + str(tot_sum_neg[key]) + '\t' + str(tot_tot_counts[key]) + '\n')

        final_counts.close()

        # Sort the file numerically
        os.system('sort -k1,2 -n ' + prefix + 'SNP_COUNTS.txt ' + ' -o ' + prefix + 'SNP_COUNTS.txt')

        # Clean up intermediate files.
        if args.noclean is False:
            os.system('rm *.err *.out *COUNTS_* *split_sam_* *.qsub')
            if wasbam is True:
                os.system('rm ' + sam_file)

    ###############
    # SINGLE MODE #
    ###############

    # If we're running in single mode (each job submitted by multiplex mode will be running in single mode)
    elif args.mode == 'single':

        # First read in the information on the SNPs that we're interested in.
        snps = {}    # Initialize a dictionary of SNP positions

        snp_file = run.open_zipped(args.snps, 'r')
        for line in snp_file:
            line = line.rstrip('\n')
            line_t = line.split('\t')

            pos = str(line_t[0]) + '|' + str(line_t[2])
            snps[pos] = line_t[3]

        snp_file.close()

        potsnp_dict = {}    # This is the dictionary of potential SNPs for each read.

        # Now parse the SAM file to extract only reads overlapping SNPs.
        in_sam = open(sam_file, 'r')
        for line in in_sam:
            if re.match('^@', line):    # Write header lines if applicable
                continue

            # Skip lines that overlap indels OR don't match Ns
            line = line.rstrip('\n')
            line_t = line.split('\t')

            if 'D' in line_t[5] or 'I' in line_t[5]:
                continue

            # Split the tags to find the MD tag:
            tags = line_t[11].split(' ')
            for i in tags:
                if re.match('^MD:', i) and 'N' in i:

                    # Remember that, for now, we're not allowing reads that overlap insertions/deletions.

                    chr = line_t[2]
                    pos = int(line_t[3])-1
                    read = line_t[9]

                    read_seq = ''

                    # Need to determine whether it's forward or reverse complimented based on the bitwise
                    # flag. This is based on the orientation bit '0b10000'  0 = forward, 1 = reverse, and
                    # the mate pairing bits (0b1000000, first mate; 0b10000000, second mate). We're assuming
                    # correct mapping such that FIRST MATES on the NEGATIVE STRAND are NEGATIVE, while
                    # SECOND MATES on the NEGATIVE STRAND are POSITIVE.

                    if args.single is True:
                        flag = int(line_t[1])
                        if flag & 0b10000:  # RYO UPDATED HERE
                            orientation = '-'
                        else:
                            orientation = '+'

                    else:
                        flag = int(line_t[1])
                        if flag & 0b1000000:    # First mate
                            if flag & 0b10000:    # If reverse, then negative strand
                                orientation = '-'
                            else:
                                orientation = '+'

                        elif flag & 0b10000000:    # Second mate
                            if flag & 0b10000:    # If reverse, then positive
                                orientation = '+'
                            else:
                                orientation = '-'

                    # Parse the CIGAR string
                    cigar_types, cigar_vals = split_CIGAR(line_t[5])

                    if cigar_types[0] == 'S':
                        MD_start = int(cigar_vals[0])
                    else:
                        MD_start = 0

                    # Get the genomic positions corresponding to each base-pair of the read
                    read_genomic_positions = CIGAR_to_Genomic_Positions(cigar_types, cigar_vals, line_t[3])

                    # Get the tag data
                    MD_vals = i.split(':')
                    MD_split = re.findall('\d+|\D+', MD_vals[2])

                    genome_start = 0

                    # The snp_pos dictionary will store the 1-base position => allele
                    snp_pos = {}
                    for i in MD_split:
                        if re.match('\^', i):
                            pass
                        elif i.isalpha():
                            if i == 'N':
                                snp_pos[read_genomic_positions[genome_start]] = read[MD_start]
                                MD_start += 1
                                genome_start += 1
                            else:
                                MD_start += 1
                                genome_start += 1
                        else:
                            MD_start += int(i)
                            genome_start += int(i)

                    for i in snp_pos:

                        # RYO: START EDIT - Implemented Filter
                        posVal = str(line_t[2]) + '|' + str(i)
                        if posVal not in snps:
                            continue
                        # RYO: END EDIT - Implmented Filter

                        snp = str(line_t[2]) + '|' + str(i) + '\t' + str(snp_pos[i]) + '\t' + orientation
                        if str(line_t[0]) in potsnp_dict:
                            if snp not in potsnp_dict[str(line_t[0])]:
                                potsnp_dict[str(line_t[0])].append(snp)  # RYO EDIT HERE - added conditional so that pairs of reads are not considered twice if they both overlap the same snp.
                        else:
                            potsnp_dict[str(line_t[0])] = []
                            potsnp_dict[str(line_t[0])].append(snp)

        in_sam.close()

        # Initialize the counting dictionaries
        pos_counts = {}
        neg_counts = {}

        # Go through the potential SNP dictionary and choose one SNP at random for those overlapping multiple SNPs
        keys = list(potsnp_dict.keys())
        for key in keys:
            snp = random.choice(potsnp_dict[key]).split('\t')

            if snp[0] in snps:
                if snp[0] in pos_counts or snp[0] in neg_counts:
                    if snp[2] == '+':
                        if snp[1] == 'A':
                            pos_counts[snp[0]][0] += 1
                        if snp[1] == 'C':
                            pos_counts[snp[0]][1] += 1
                        if snp[1] == 'G':
                            pos_counts[snp[0]][2] += 1
                        if snp[1] == 'T':
                            pos_counts[snp[0]][3] += 1

                    elif snp[2] == '-':
                        if snp[1] == 'A':
                            neg_counts[snp[0]][0] += 1
                        if snp[1] == 'C':
                            neg_counts[snp[0]][1] += 1
                        if snp[1] == 'G':
                            neg_counts[snp[0]][2] += 1
                        if snp[1] == 'T':
                            neg_counts[snp[0]][3] += 1

                else:
                    pos_counts[snp[0]] = [0, 0, 0, 0]
                    neg_counts[snp[0]] = [0, 0, 0, 0]
                    if snp[2] == '+':
                        if snp[1] == 'A':
                            pos_counts[snp[0]][0] += 1
                        if snp[1] == 'C':
                            pos_counts[snp[0]][1] += 1
                        if snp[1] == 'G':
                            pos_counts[snp[0]][2] += 1
                        if snp[1] == 'T':
                            pos_counts[snp[0]][3] += 1

                    elif snp[2] == '-':
                        if snp[1] == 'A':
                            neg_counts[snp[0]][0] += 1
                        if snp[1] == 'C':
                            neg_counts[snp[0]][1] += 1
                        if snp[1] == 'G':
                            neg_counts[snp[0]][2] += 1
                        if snp[1] == 'T':
                            neg_counts[snp[0]][3] += 1

        # Open the output file and write the SNP counts to it

        if not args.suffix:
            out_counts = run.open_zipped(prefix + 'SNP_COUNTS.txt', 'w')
        else:
            out_counts = run.open_zipped(prefix + 'SNP_COUNTS_' + args.suffix, 'w')

        # Write header
        out_counts.write('CHR\tPOSITION\tPOS_A|C|G|T\tNEG_A|C|G|T\tSUM_POS_READS\tSUM_NEG_READS\tSUM_READS\n')

        # Sort SNP positions and write them
        keys = list(pos_counts.keys())
        keys.sort()

        for key in keys:
            pos = key.split('|')
            sum_pos = sum(pos_counts[key])
            sum_neg = sum(neg_counts[key])
            tot_sum = sum(pos_counts[key]) + sum(neg_counts[key])
            pos_fix = [str(x) for x in pos_counts[key]]
            neg_fix = [str(x) for x in neg_counts[key]]
            positive = '|'.join(pos_fix)
            negative = '|'.join(neg_fix)

            out_counts.write(pos[0] + '\t' + pos[1] + '\t' + positive + '\t' + negative + '\t' + str(sum_pos) + '\t' + str(sum_neg) + '\t' + str(tot_sum) + '\n')

        out_counts.close()

        if not args.suffix:
            pass
        else:
            os.system('touch ' + prefix + args.suffix + '_done')

if __name__ == '__main__' and '__file__' in globals():
    sys.exit(main())
