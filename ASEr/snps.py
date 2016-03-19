"""
Filter SNP lists by individual.

============================================================================

        AUTHOR: Michael D Dacre, mike.dacre@gmail.com
  ORGANIZATION: Stanford University
       LICENSE: MIT License, property of Stanford, use as you wish
       VERSION: 0.1
       CREATED: 2016-36-14 12:03
 Last modified: 2016-03-18 17:05

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
#                               Useful Classes                                #
###############################################################################


class Individual(object):

    """A simple individual to hold a frozenset of snps."""

    def add_bed(self, bedfile):
        """Add a list of pybedtools Interval objects to self as self.bed.

        Requires pybedtools, adds only records for snps in this individual.

        Note: This is a slow operation.
        """
        try:
            from pybedtools import BedTool
        except ImportError:
            logme.log('pybedtools is not installed.\n' +
                      'Please install and try again. You can get it from here:\n' +
                      'https://github.com/daler/pybedtools',
                      level='error')
            return -1
        bed = BedTool(bedfile)
        self.bed = [i for i in bed.filter(lambda a: a.name in self.snps)]

    def save_bed(self, outfile, bedfile=None):
        """Save a bed file of SNPs to outfile.

        :bedfile: Only required if self.bed has not been populated already.
        """
        if not self.bed:
            if bedfile:
                self.add_bed(bedfile)
            else:
                logme.log('Cannot save a bed file with no starting bed.\n' +
                          'Try again with the bedfile option, or by\n' +
                          'running add_bed first.', level='error')
                return -1
        with open_zipped(outfile, 'w') as fout:
            for i in self.bed:
                fout.write(str(i))
        return 0

    def __init__(self, name, snplist):
        """Store the name and snplist."""
        self.name = name
        self.snps = frozenset(snplist)

    def __getattr__(self, attr):
        """Set length if not already set."""
        if attr == '_len':
            self._len = len(self.snps)
            return self._len

    def __repr__(self):
        """Print information about the individual."""
        return '<Individual(Name={},SNPs={})>'.format(self.name, self._len)

    def __str__(self):
        """Print name."""
        return self.name

    def __len__(self):
        """Length of snp list."""
        return self._len

    def __contains__(self, item):
        """Check for item in snps."""
        return item in self.snps

    def __iter__(self):
        """Iterate through snps."""
        for snp in self.snps:
            yield snp

###############################################################################
#                                File Filters                                 #
###############################################################################


def filter_bed(bedfile, snp_list, outfile=sys.stdout):
    """Filter a bedfile to only include snps in snp_list, print to outfile.

    :bedfile:  A bed file of all the SNPs, can be gzipped.
    :snp_list: List/tuple/set/frozenset of snp names.
    :outfile:  Something .bed or .bed.gz, deault STDOUT.
    :returns:  0 on success 1 on failure

    """
    try:
        from pybedtools import BedTool
    except ImportError:
        logme.log('pybedtools is not installed.\n' +
                  'Please install and try again. You can get it from here:\n' +
                  'https://github.com/daler/pybedtools',
                  level='error')
        return -1

    if not isinstance(snp_list, (tuple, list, set, frozenset)):
        raise Exception('snp_list must be tuple/list/set/frozenset ' +
                        'it is: {}'.format(type(snp_list)))

    bed      = BedTool(bedfile)
    filtered = bed.filter(lambda a: a.name in snp_list)

    with open_zipped(outfile, 'w') as fout:
        fout.write(str(filtered))


def get_het_snps_from_recodeAD(infile, snps=None, individuals=None,
                               split_individual=None, name_index=0):
    """Iterator to return a list of SNPs one individual at a time.

    Yields: An Individual class with .name, and .snps.

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
        if not isinstance(snps, frozenset):
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
                shortname = fields[0].split(split_individual)[name_index] \
                    if split_individual else fields[0]
                if shortname not in individuals:
                    continue
            else:
                shortname = name
            hets = frozenset([i for i, j in zip(headers, fields[7::2])
                              if j == '1'])
            if snps:
                hets = hets.intersection(snps)

            yield Individual(shortname, hets)


def snps_from_bed(snp_file):
    """Return a frozenset of SNP names from a bed file."""
    snps = []
    with open_zipped(snp_file) as fin:
        for i in fin:
            snps.append(i.split('\t')[3])
    return frozenset(snps)


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
