#!/usr/bin/env python
# -*- coding: utf-8 -*-


# Created on: 2015.03.16
# Author: Carlo Artieri

###########
# MODULES #
###########
import sys              # Access to simple command-line arguments
import argparse         # Access to long command-line parsing
import datetime         # Access to calendar/clock functions
import re               # Access to REGEX splitting
import math             # Access to math functions
import random           # Access to random number generation
import subprocess       # Access to direct command line in/out processing
import os
import textwrap         # Add text block wrapping properties
from time import sleep  # Allow system pausing

# US
from ASEr.run import open_zipped

##########################
# COMMAND-LINE ARGUMENTS #
##########################

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


class CustomFormatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawDescriptionHelpFormatter):
    pass

parser = argparse.ArgumentParser(description='This script will take a BAM file mapped to a SNP-masked genome and count the number of reads overlapping each SNP. When reads overlap multiple SNPs, the read is randomly assigned to one of them.', add_help=False, epilog=epilog, formatter_class=CustomFormatter)
req = parser.add_argument_group('Required arguments')
req.add_argument('-m', '--mode', action='store', dest='mode', help='Operation mode', choices=['single', 'multi'], required=True, metavar='mode')
req.add_argument('-s', '--snps', action='store', dest='snps', help='SNP BED file', required=True, metavar='<BED>')
req.add_argument('-r', '--reads', action='store', dest='reads', help='Mapped reads file [sam or bam]', required=True, metavar='<[S/B]AM>')
uni = parser.add_argument_group('Universal optional arguments')
uni.add_argument('-p', '--prefix', action='store', dest='prefix', help='Prefix for temp files and output', default='TEST', metavar='')
uni.add_argument('-b', '--bam', action='store_true', dest='bam', help='Mapped read file type is bam (auto-detected if *.bam)')
uni.add_argument('-t', '--single', action='store_true', dest='single', help='Mapped reads are single-end')
uni.add_argument('-n', '--noclean', action='store_true', dest='noclean', help='Do not delete intermediate files (for debuging)')
uni.add_argument('-h', '--help', action='help', help='show this help message and exit')
mult = parser.add_argument_group('Multi(plex) mode arguments')
mult.add_argument('-j', '--jobs', action='store', dest='jobs', type=int, help='Divide into # of jobs', default=100, metavar='')
mult.add_argument('-w', '--walltime', action='store', dest='walltime', help='Walltime for each job', default='3:00:00', metavar='')
mult.add_argument('-k', '--mem', action='store', dest='memory', help='Memory for each job', default='5000MB', metavar='')
single = parser.add_argument_group('Single mode arguments')
single.add_argument('-f', '--suffix', action='store', dest='suffix', help='Suffix for multiplexing [set automatically]', default='', metavar='')

args = parser.parse_args()

#############
# FUNCTIONS #
#############


# Convert a FASTA file to a dictionary where keys = headers and values are the sequence
def fasta_to_dict(file):

    fasta_file = open_zipped(file, "r")    # Open the file for reading
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


# Split the CIGAR string and return two lists: types and values in order
def split_CIGAR(cigar):

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


# Using the CIGAR string return a list of genomic coordinates corresponding to the
# individual bases of the read to get SNP positions from the MD tag.
def CIGAR_to_Genomic_Positions(cigar_types, cigar_vals, pos):

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
            genomic_positions = genomic_positions + list(range(int(curr_pos), int(curr_pos)+int(cigar_vals[i])))
            curr_pos = int(curr_pos) + int(cigar_vals[i])
    return genomic_positions

##########
# SCRIPT #
##########

# Initialize variables
prefix = args.prefix + '_'
wasbam = False

# Check if the read file is sam or bam
file_check = args.reads.split('.')
file_check[-1] = file_check[-1].lower()
if file_check[-1] == 'bam' or args.bam is True:
    wasbam = True
    sam_path, sam_file = os.path.split('.'.join(args.reads.split('.')[:-1]) +
                                       '.sam')
    sys.stderr.write('Converting BAM to SAM ...\n')
    os.system('samtools view ' + args.reads + ' > ' + sam_file)
    sys.stderr.write('Done\n')

elif file_check[-1] == 'sam' or args.bam is False:
    sam_file = args.reads


##################
# MULTIPLEX MODE #
##################

# If we're running in multiplex mode
if args.mode == 'multi':

    # Determine how many reads will be in each split sam file.
    num_lines = os.popen('wc -l ' + sam_file + ' | awk \'{print $1}\'').read()
    num_reads = int(int(num_lines)/args.jobs) + 1

    # Subset the SAM file into X number of jobs
    cnt = 0
    currjob = 1
    suffix = '.split_sam_' + str(currjob).zfill(4)
    run_file = os.path.join(sam_path, prefix + sam_file + suffix)

    sam_split = open(run_file, 'w')

    in_sam = open(sam_file, 'r')
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
                sam_split = open(prefix + sam_file + suffix, 'w')
                sam_split.write(line2)
                cnt = 0

    in_sam.close()
    sam_split.close()

    # Create PBS scripts and submit jobs to the cluster

    if args.noclean is True:
        subnoclean = '--noclean'
    else:
        subnoclean = ''

    for i in range(1, args.jobs+1):
        suffix = str(i).zfill(4)
        reads_file = prefix + sam_file + '.split_sam_' + suffix

        # qsub script modify as necessary
        qsub_script = """\
        # PBS -m n
        # PBS -V
        # PBS -d ./
        # PBS -N """ + prefix + suffix + """
        # PBS -l nodes=1:ppn=1
        # PBS -l walltime=""" + args.walltime + """
        # PBS -l mem=""" + args.memory + """
        # PBS -e """ + prefix + suffix + """_err.txt
        # PBS -o """ + prefix + suffix + """_out.txt
        python2 """ + parser.prog + """ --mode single --snps """ + args.snps + """ --reads """ + \
            reads_file + """ --suffix """ + suffix + """ --prefix """ + args.prefix + """
        exit 0
        """
        qsub = open('qsub.txt', 'w')
        qsub.write(textwrap.dedent(qsub_script))
        qsub.close()

        # Submit jobs to queue
        os.system('qsub qsub.txt')
        sleep(2)    # Pause for two seconds to make sure job is properly submitted

    # Now wait and check for all jobs to complete every so long
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

    # Once the jobs are done, concatenate all of the counts into one file.
    # Initialize dictionaries

    os.system('rm *_done')    # Remove the 'done' files in case we want to run again.

    tot_pos_counts = {}
    tot_neg_counts = {}
    tot_tot_counts = {}
    tot_sum_pos = {}
    tot_sum_neg = {}

    for i in range(1, args.jobs+1):
        suffix = str(i).zfill(4)
        in_counts = open_zipped(prefix + 'SNP_COUNTS_' + suffix, 'r')

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
    final_counts = open_zipped(prefix + 'SNP_COUNTS.txt', 'w')
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
        os.system('rm *err.txt *out.txt *COUNTS_* *split_sam_* qsub.txt')
        if wasbam is True:
            os.system('rm ' + sam_file)

###############
# SINGLE MODE #
###############

# If we're running in single mode (each job submitted by multiplex mode will be running in single mode)
elif args.mode == 'single':

    # First read in the information on the SNPs that we're interested in.
    snps = {}    # Initialize a dictionary of SNP positions

    snp_file = open_zipped(args.snps, 'r')
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
        out_counts = open_zipped(prefix + 'SNP_COUNTS.txt', 'w')
    else:
        out_counts = open_zipped(prefix + 'SNP_COUNTS_' + args.suffix, 'w')

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
