"""
Filter SNP lists by individual.

============================================================================

        AUTHOR: Michael D Dacre, mike.dacre@gmail.com
  ORGANIZATION: Stanford University
       LICENSE: MIT License, property of Stanford, use as you wish
       VERSION: 0.1
       CREATED: 2016-36-14 12:03
 Last modified: 2016-03-17 12:08

   DESCRIPTION:

============================================================================
"""
import sys
from subprocess import check_output

# Our functions
from .plink import is_recodeAD
from .plink import PlinkError
from .run import open_zipped
from .run import is_file_type
from .run import write_iterable

# Logging
from . import logme
logme.LOGFILE   = sys.stderr
logme.MIN_LEVEL = 'info'

__all__ = ['get_het_snps_from_recodeAD', 'filter_snps_by_exon']


###############################################################################
#                                File Filters                                 #
###############################################################################


def get_het_snps_from_recodeAD(infile, snps=None, individuals=None,
                               split_individual=None, name_index=0):
    """Iterator to return a list of SNPs one individual at a time.

    Requires a raw recodeAD file from plink.

    Returns a set of all heterozygous SNPs for each individual.

    :infile:           The recodeAD file
    :snps:             An list, set, or tuple of SNPs to filter against,
                       only SNPs that are in this list will be returned.
    :individuals:      A list, set, or tuple of individuals.
    :split_individual: Split the individual name by this character.
    :name_index:       If split_individual used, use this index to choose the
                       name element.
    """
    if not is_recodeAD(infile):
        raise PlinkError('{} is not a recodeAD file'.format(infile))

    # Check data types
    if snps:
        if isinstance(snps, (list, tuple, set)):
            snps = frozenset(snps)
        elif isinstance(snps, frozenset):
            snps = snps
        else:
            logme.log('snps is wrong type', 'critical')
            raise TypeError('snps is {}, '.format(type(snps)) +
                            'must be list, tuple, set, or frozenset')
    if individuals:
        if isinstance(individuals, (list, tuple, set)):
            individuals = frozenset(individuals)
        elif isinstance(individuals, frozenset):
            individuals = individuals
        else:
            logme.log('individuals is wrong type', 'critical')
            raise TypeError('individuals is {}, '.format(type(individuals)) +
                            'must be list, tuple, set, or frozenset')
    if split_individual and len(split_individual) > 1:
        raise Exception('split_individual must be a single character')
    name_index = int(name_index)

    # Parse the file
    with open_zipped(infile) as fin:
        headers = [i[:-4] for i in fin.readline().rstrip().split(' ')[7::2]]
        for line in fin:
            fields = line.rstrip().split(' ')
            name = fields[0]
            if individuals:
                name2 = fields[0].split(split_individual)[name_index] \
                    if split_individual else fields[0]
                if name2 not in individuals:
                    continue
            hets = frozenset([i for i, j in zip(headers, fields[7::2])
                              if j == '1'])
            if snps:
                yield hets.intersection(snps)
            else:
                yield hets


def filter_snps_by_exon(snp_file, exon_file, outfile=None, outbed=False):
    """Filter a bed file of SNPs and return frozenset of SNPs in exons.

    NOTE: This module uses pybedtools and will abort if it isn't installed.

    :snp_file:  bed/gff/gtf file of SNPs (gzipped is OK)
    :exon_file: bed/gff/gtf file of exons (gzipped is OK)
    :outfile:   If provided, the resultant SNPs are saved to that file, if
                file ends in bed the reads will be saved in that
                format, any other ending will be saved as a list of reads.
                Note: .gz does not count as an ending and just defines
                compression (i.e. bed.gz will be a gzipped bed file)
    :outbed:    If provided, bedfile output is forced
    """
    try:
        from pybedtools import BedTool
    except ImportError:
        logme.log('pybedtools is not installed.\n' +
                  'Please install and try again. You can get it from here:\n' +
                  'https://github.com/daler/pybedtools',
                  level='error')
        return -1

    snps  = BedTool(snp_file)
    exons = BedTool(exon_file)

    intersection = snps.intersect(exons)

    final_snps = frozenset([i.name for i in intersection])

    if outfile:
        if is_file_type(outfile, 'bed') or outbed:
            with open_zipped(outfile, 'w') as fout:
                fout.write(str(intersection))
        else:
            write_iterable(final_snps, outfile)

    return final_snps
