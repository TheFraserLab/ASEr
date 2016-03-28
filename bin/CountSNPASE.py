#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Count number of reads overlapping each SNP in a sam/bam file.

============================================================================

        AUTHOR: Carlo Artieri
    MAINTAINER: Michael D Dacre, mike.dacre@gmail.com
  ORGANIZATION: Stanford University
       LICENSE: MIT License, property of Stanford, use as you wish
       CREATED: 2015-03-16
 Last modified: 2016-03-24 14:45

   DESCRIPTION: This script will take a BAM file mapped to a SNP-masked
                genome and count the number of reads overlapping each SNP.
                When reads overlap multiple SNPs, the read is randomly
                assigned to one of them.

============================================================================
"""
import os                  # Path manipulation
import sys                 # Access to simple command-line arguments
import argparse            # Access to long command-line parsing
import re                  # Access to REGEX splitting
import random              # Access to random number generation
from time import sleep     # Allow system pausing
from multiprocessing import cpu_count
from pysam import Samfile  # Read sam and bamfiles

# Us
from ASEr import logme     # Logging functions
from ASEr import run       # File handling functions
from ASEr import cluster   # Queue submission
from ASEr.snps import chrom_to_num  # Chromosome number standardization

# Logging
logme.MIN_LEVEL = 'info'  # Switch to 'debug' for more verbose, 'warn' for less

###############################################################################
#                          Command Line Description                           #
###############################################################################

EPILOG = """\
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

###############################################################################
#                                  Functions                                  #
###############################################################################


###
#CARLO'S COMMENTS
#
#I can't remember when/if this function gets used, but I've since learned of a MUCH faster way
#to access a FASTA sequence without storing the whole stupid thing in memory, which will kill
#most applications on the human genome.
#
#First, the FASTA must be indexed with samtools:
#
# samtools faidx ref.fasta
#
#Then we can use pysam to return any 0-based sequence as follows:
#
# with pysam.FastaFile(ref.fasta, 'r') as in_fasta:
#     seq = in_fasta.fetch(reference=CHROM, start=START, end=END)
#
#Because it's indexed and uses a C-API, it's very fast, and low mem.

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
    """Split the CIGAR string and return two lists: types and values."""
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


def count_reads(samfilename):
    """Return the number of reads in a samfile.

    Use built-in counts if possible, but fall back to looping through otherwise
    """
    sf = Samfile(samfilename)
    try:
        return sf.mapped
    except ValueError:
        num_reads = 0
        for num_reads, _ in enumerate(sf):
            pass
    return num_reads


def split_samfile(sam_file, splits, prefix='', path=''):
    """Take a sam file and split it splits number of times.

    :path:    Where to put the split files.
    :prefix:  A prefix for the outfile names.
    :returns: A tuple of job files.
    """
    # Determine how many reads will be in each split sam file.
    num_lines = count_reads(sam_file)
    num_reads = int(int(num_lines)/splits) + 1

    # Get rid of starting path
    sam_name = os.path.basename(sam_file)

    # Subset the SAM file into X number of jobs
    cnt      = 0
    currjob  = 1
    suffix   = '.split_sam_' + str(currjob).zfill(4)
    run_file = os.path.join(path, prefix + sam_name + suffix)
    rmode    = 'rb' if sam_name.split('.')[0] == 'bam' else 'r'
    wmode    = 'wb'

    # Actually split the file
    outfiles = [run_file]
    with Samfile(sam_file, rmode) as in_sam:
        sam_split = Samfile(run_file, wmode, template=in_sam)
        for line in in_sam:
            cnt += 1
            if cnt < num_reads:
                sam_split.write(line)
            elif cnt == num_reads:
                # Check if next line is mate-pair. If so, don't split.
                line2     = next(in_sam)
                currjob  += 1
                suffix    = '.split_sam_' + str(currjob).zfill(4)
                run_file  = os.path.join(path, prefix + sam_name + suffix)
                new_sam   = Samfile(run_file, wmode, template=in_sam)
                outfiles.append(run_file)

                if line.qname == line2.qname:
                    sam_split.write(line)
                    sam_split.write(line2)
                    sam_split.close()
                    cnt = 0
                else:
                    sam_split.write(line)
                    sam_split.close()
                    new_sam.write(line2)
                    cnt = 0
                sam_split = new_sam
        sam_split.close()
    return tuple(outfiles)


###############################################################################
#                                 Main Script                                 #
###############################################################################


def main(argv=None):
    """Main script."""
    ##########################
    # COMMAND-LINE ARGUMENTS #
    ##########################

    # Get myself
    program_name = sys.argv[0]

    if not argv:
        argv = sys.argv[1:]

    # Get the cluster type used to control arguments
    cluster_type = cluster.get_cluster_environment()

    parser  = argparse.ArgumentParser(
        description=__doc__, add_help=False,
        epilog=EPILOG, formatter_class=run.CustomFormatter)

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
    uni.add_argument('-n', '--noclean', action='store_true',
                     help='Do not delete intermediate files (for debuging)')
    uni.add_argument('-R', '--random-seed', default=None, type=int,
                     help='Set the state of the randomizer (for testing)')
    uni.add_argument('-h', '--help', action='help',
                     help='show this help message and exit')

    mult = parser.add_argument_group('Multi(plex) mode arguments')
    mult.add_argument('-j', '--jobs', type=int,
                      help='Divide into # of jobs', default=100, metavar='')
    if cluster_type == 'slurm' or cluster_type == 'torque':
        mult.add_argument('-w', '--walltime',
                          help='Walltime for each job', default='3:00:00',
                          metavar='')
        mult.add_argument('-k', '--mem', dest='memory', metavar='',
                          help='Memory for each job', default='5000MB')
        mult.add_argument('--queue',
                          help='Queue to submit jobs to', default='batch',
                          metavar='')
        mult.add_argument('--cluster', choices=['torque', 'slurm', 'normal'],
                          help='Which cluster to use, normal uses threads ' +
                          'on this machine', default=cluster_type)
    mult.add_argument('--threads', type=int, metavar='', default=cpu_count(),
                      help='Max number of threads to run at a time ' +
                      '(normal mode only).')

    single = parser.add_argument_group('Single mode arguments')
    single.add_argument('-f', '--suffix', default='', metavar='',
                        help='Suffix for multiplexing [set automatically]')

    logging = parser.add_argument_group('Logging options')
    logging.add_argument('-q', '--quiet', action='store_true',
                         help="Quiet mode, only prints warnings.")
    logging.add_argument('-v', '--verbose', action='store_true',
                         help="Verbose mode, prints debug info too.")
    logging.add_argument('--logfile',
                         help='Logfile to write messages too, default is ' +
                         'STDERR')

    args = parser.parse_args()
    if args.random_seed is not None:
        random.seed(args.random_seed)
        print("Seed: ", args.random_seed, random.getstate()[1][:10])

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

    if args.reads.endswith('bam') or args.bam:
        mode = 'rb'
    else:
        mode = 'r'

    ##################
    # MULTIPLEX MODE #
    ##################

    # If we're running in multiplex mode
    if args.mode == 'multi':
        logme.log('Splitting sam file {} into {} files.'.format(sam_file,
                                                                args.jobs))
        reads_files = split_samfile(os.path.join(sam_path, sam_file),
                                    args.jobs, prefix)
        logme.log('Splitting complete.')

        # Create PBS scripts and submit jobs to the cluster
        subnoclean = ' --noclean' if args.noclean else ''
        logme.log('Submitting split files to cluster')
        jobs = []  # Hold job info for later checking
        for reads_file in reads_files:
            suffix = reads_file[-4:]

            command = ("python2 " + program_name + " --mode single --snps " +
                       args.snps + " --reads " + reads_file + " --suffix " +
                       suffix + " --prefix " + args.prefix + subnoclean +
                       ' --bam')

            if cluster_type == 'normal':
                jobs.append(cluster.submit(command, name=prefix + suffix,
                                           threads=args.threads))
            else:
                jobs.append(cluster.submit(command, name=prefix + suffix,
                                           time=args.walltime, cores=1,
                                           mem=args.memory,
                                           partition=args.queue))
            sleep(2)    # Pause for two seconds to make sure job is submitted

        # Now wait and check for all jobs to complete every so long
        logme.log('Submission done, waiting for jobs to complete.')

        # First wait for jobs in queue to complete
        cluster.wait(jobs)
        sleep(1)

        # Next, check if any jobs failed
        failed = []
        for i in range(1, args.jobs+1):
            suffix = str(i).zfill(4)
            if not os.path.isfile(prefix + suffix + '_done'):
                failed.append(prefix + suffix)

        # If any jobs failed, terminate
        if failed:
            logme.log('Some jobs failed!', 'critical')
            return -1

        logme.log('Jobs completed.')
        # Remove 'done' files in case we want to run again.
        os.system('rm {prefix}*_done'.format(prefix=prefix))

        # Once the jobs are done, concatenate all of the counts into one file.
        # Initialize dictionaries

        tot_pos_counts = {}
        tot_neg_counts = {}
        tot_tot_counts = {}
        tot_sum_pos = {}
        tot_sum_neg = {}

        for i in range(1, args.jobs+1):
            suffix = str(i).zfill(4)
            in_counts = prefix + 'SNP_COUNTS_' + suffix

            # Parse the line to add it to the total file
            with run.open_zipped(in_counts, 'r') as in_counts:
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

        # Write out the final concatenated file
        with run.open_zipped(prefix + 'SNP_COUNTS.txt', 'w') as final_counts:
            final_counts.write('CHR\tPOSITION\tPOS_A|C|G|T\tNEG_A|C|G|T\t' +
                               'SUM_POS_READS\tSUM_NEG_READS\tSUM_READS\n')

            keys = sorted(tot_pos_counts.keys())

            for key in keys:
                pos = key.split('|')
                pos_fix = [str(x) for x in tot_pos_counts[key]]
                neg_fix = [str(x) for x in tot_neg_counts[key]]
                pos_out = '|'.join(pos_fix)
                neg_out = '|'.join(neg_fix)
                final_counts.write(str(pos[0]) + '\t' + str(pos[1]) + '\t' +
                                   pos_out + '\t' + neg_out + '\t' +
                                   str(tot_sum_pos[key]) + '\t' +
                                   str(tot_sum_neg[key]) + '\t' +
                                   str(tot_tot_counts[key]) + '\n')

        # Sort the file numerically
        os.system('sort -k1,2 -n ' + prefix + 'SNP_COUNTS.txt ' + ' -o ' +
                  prefix + 'SNP_COUNTS.txt')

        # Clean up intermediate files.
        if args.noclean is False:
            cluster.clean()
            os.system('rm {prefix}*COUNTS_* {prefix}*split_sam_*'.format(
                prefix=prefix))

    ###############
    # SINGLE MODE #
    ###############

    # If we're running in single mode (each job submitted by multiplex mode
    # will be running in single mode)
    elif args.mode == 'single':

        # First read in the information on the SNPs that we're interested in.
        snps = {}    # Initialize a dictionary of SNP positions

        with run.open_zipped(args.snps) as snp_file:
            for line in snp_file:
                line = line.rstrip('\n')
                line_t = line.split('\t')

                pos = chrom_to_num(line_t[0]) + '|' + str(line_t[2])
                snps[pos] = line_t[3]

        ###
        #CARLO'S COMMENTS
        #
        #Reading in BED/GTF/VCF files is a situation where pandas can make code easier to 
        #maintain/understand. It's not necessarily fewer lines, but it's easier to see what's
        #going on. Also, when I wrote this script, I didn't realize that tuples can act as 
        #dictionary indices, which is MUCH better than the CHROM + '|' + POS construction.
        #
        #From working with computer programmers, I've learned that, as a general rule, it's bad
        #form to use elements from split strings without first assigning them to variables. 
        #This is because it requires everyone to know the exact formating that args.snps will 
        #have. If they're first assigned to variables, or have columns indicated with pandas,
        #as below, it's much clearer to everyone what's happening.
        #
        #(assuming import pandas as pd) 
        #
        # snps = {}
        # snp_file = pd.read_table(args.bed)
        # snp_file.columns = (["CHROM","START","END","SNP"])
        # for idx,row in snp_file.itterows():
        #     snps[(chrom_to_num(row["CHROM"]),int(row["END"]))] = row["SNP"] 
        #
        ###

        # This is the dictionary of potential SNPs for each read.
        potsnp_dict = {}

        # Now parse the SAM file to extract only reads overlapping SNPs.
        in_sam     = Samfile(args.reads, mode)
        references = in_sam.references  # Faster to make a copy of references.

        # Trackers to count how many reads are lost at each step
        indel_skip = 0
        nosnp_skip = 0
        count      = 0
        snp_count  = 0
        ryo_filter = 0

        ###
        #CARLO'S COMMENTS
        #
        #This could be:
        #
        # with pysam.Samfile(args.reads,mode) as in_sam:
        #     for line in in_sam.fetch(until_eof=True): 

        for line in in_sam: 
            count += 1

            # Skip lines that overlap indels OR don't match Ns
            cigarstring = line.cigarstring

            if 'D' in cigarstring or 'I' in cigarstring:
                indel_skip += 1
                continue

            ###
            #CARLO'S COMMENTS
            #
            #The newer version of pysam allows you to get variants from reads without having
            #to manually parse tags and would replace a large amount of the code below
            #
            # for (read_pos,genome_pos,genome_seq) in read.get_aligned_pairs(with_seq=True):
            #     if genome_pos: #Is this base mapped?
            #         if ((line.reference_name,genome_pos)) in snps: #Assuming we've converted to tuples
            #
            #             #This would be a good point to include a minimum quality score filter.
            #             if line.query_quality[read_pos] >= MIN_QUALITY:
            #
            #                 #The read sequence of the SNP will be:
            #                 line.query_sequence[read_pos]
            #                 #and can be added to the appropriate counter.
            ###

            # Split the tags to find the MD tag:
            tags = line.tags
            for tagname, tagval in tags:
                if tagname == 'MD' and 'N' in tagval:
                    # Remember that, for now, we're not allowing reads that
                    # overlap insertions/deletions.

                    chrom = references[line.rname]
                    pos   = line.pos
                    read  = line.seq

                    # We're assuming
                    # correct mapping such that FIRST MATES on the NEGATIVE
                    # STRAND are NEGATIVE, while SECOND MATES on the NEGATIVE
                    # STRAND are POSITIVE.

                    if line.is_reverse:
                        orientation = '-'
                    else:
                        orientation = '+'

                    # Parse the CIGAR string
                    cigar_types, cigar_vals = split_CIGAR(cigarstring)

                    if cigar_types[0] == 'S':
                        MD_start = int(cigar_vals[0])
                    else:
                        MD_start = 0

                    # Get the genomic positions corresponding to each base-pair
                    # of the read
                    read_genomic_positions = CIGAR_to_Genomic_Positions(
                        cigar_types, cigar_vals, line.pos+1)

                    # Get the tag data
                    MD_split = re.findall('\d+|\D+', tagval)

                    genome_start = 0

                    # The snp_pos dictionary will store the 1-base position
                    # => allele
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
                        snp_count += 1

                        # RYO: START EDIT - Implemented Filter
                        posVal = line.reference_name + '|' + str(i)
                        if posVal not in snps:
                            nosnp_skip += 1
                            continue
                        # RYO: END EDIT - Implmented Filter

                        snp = '{chr}|{i}\t{snp_pos}\t{orientation}'.format(
                            chr=chrom, i=i, snp_pos=snp_pos[i],
                            orientation=orientation)
                        if line.qname in potsnp_dict:
                            if snp not in potsnp_dict[line.qname]:
                                # RYO EDIT HERE - added conditional so that
                                # pairs of reads are not considered twice if
                                # they both overlap the same snp.
                                potsnp_dict[line.qname].append(snp)
                            else:
                                ryo_filter += 1
                        else:
                            potsnp_dict[line.qname] = []
                            potsnp_dict[line.qname].append(snp)

        in_sam.close()

        # Log all of the skipped reads
        logme.log('Total reads: {}'.format(count), 'debug')
        logme.log('Reads skipped for indels: {}'.format(indel_skip), 'debug')
        logme.log('Total SNPs checked: {}'.format(snp_count), 'debug')
        logme.log('SNPs not in SNP list: {}'.format(nosnp_skip), 'debug')
        logme.log('Ryo filter: {}'.format(ryo_filter), 'debug')

        # Initialize the counting dictionaries
        pos_counts = {}
        neg_counts = {}

        # Go through the potential SNP dictionary and choose one SNP at random
        # for those overlapping multiple SNPs
        if args.random_seed is not None:  # Dictionaries are unordered, so must sort for consistent random seed output.
            keys = sorted(list(potsnp_dict.keys()))
        else:  # Because sorting is slow, only do it if random seed is set, slowdown is about 0.1s per 1 million reads..
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

        out_counts = prefix + 'SNP_COUNTS_' + args.suffix if args.suffix \
            else prefix + 'SNP_COUNTS.txt'

        with open(out_counts, 'w') as out_counts:
            # Write header
            out_counts.write('CHR\tPOSITION\tPOS_A|C|G|T\tNEG_A|C|G|T\t' +
                             'SUM_POS_READS\tSUM_NEG_READS\tSUM_READS\n')

            # Sort SNP positions and write them
            keys = sorted(pos_counts.keys())

            for key in keys:
                pos = key.split('|')
                sum_pos = sum(pos_counts[key])
                sum_neg = sum(neg_counts[key])
                tot_sum = sum(pos_counts[key]) + sum(neg_counts[key])
                pos_fix = [str(x) for x in pos_counts[key]]
                neg_fix = [str(x) for x in neg_counts[key]]
                positive = '|'.join(pos_fix)
                negative = '|'.join(neg_fix)

                out_counts.write(pos[0] + '\t' + pos[1] + '\t' + positive +
                                 '\t' + negative + '\t' + str(sum_pos) + '\t' +
                                 str(sum_neg) + '\t' + str(tot_sum) + '\n')

        if args.suffix:
            os.system('touch ' + prefix + args.suffix + '_done')

if __name__ == '__main__' and '__file__' in globals():
    sys.exit(main())
