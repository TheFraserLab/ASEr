#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Test the ASEr pipeline for sanity.

============================================================================

        AUTHOR: Michael D Dacre, mike.dacre@gmail.com
  ORGANIZATION: Stanford University
       LICENSE: MIT License, property of Stanford, use as you wish
       CREATED: 2016-48-24 10:03
 Last modified: 2016-03-24 12:58

   DESCRIPTION: The logic of this is:
                    - Download known reference samples
                    - Run the pipeline with a set random seed (0)
                    - Compare the output to a hash of a known file
                Test are passed if hash tests are all OK, and all supporting
                tests run fine.

         USAGE: This file is intended to be run with pytest from the root
                directory. It should be compatible with travis-ci

============================================================================
"""
import os
import sys
import hashlib

# Python 2/3 compatibility
try:
    from urllib2 import urlopen
except ImportError:
    from urllib.request import urlopen

# Us
from ASEr import run


#################
#  File Hashes  #
#################


# Using SHA1 hashes, SHORT hashes only include first 65,000 bytes of file.
ROOT_URL = 'http://web.stanford.edu/~pcombs/asetest/'
FILES    = {'ase.bam':        {'url': ROOT_URL + 'ase.bam',
                               'short_hash':
                               '507daac84c256c4188226366bb944e87858e0adb'},
            'variants.bed':   {'url': ROOT_URL + 'variants.bed',
                               'short_hash':
                               'cb2d7c551465f176ee3c8b4d174bdaac199a8c31'},
            'ref.gtf':        {'url': ROOT_URL + 'ref.gtf',
                               'short_hash':
                               'b9f65d2c9a274112829590bebc0b4b6de9560dc1'},
            'SNP_COUNTS.txt': {'url': ROOT_URL + 'reference_SNP_COUNTS.txt',
                               'hash':
                               '3e583f4d45e177e88fd70c4099476e9f3572fe5e'},
            'gene_ase.tsv':   {'url': ROOT_URL + 'reference_gene_ase.tsv',
                               'hash':
                               'a3589a8653c22246a07dd0fc4c582f5ad4fd493e'}
            }
ROOT_DIR = '.'
TEST_DIR = 'testdir'

###############################################################################
#                              General Functions                              #
###############################################################################

#############
#  Hashing  #
#############


def hash_file(infile):
    """Return a hash of the infile."""
    hasher = hashlib.sha1()
    with open(infile, 'rb') as fin:
        hasher.update(fin.read())
    return str(hasher.hexdigest())


def quick_hash(infile):
    """Return a hash of the infile."""
    hasher = hashlib.sha1()
    with open(infile, 'rb') as fin:
        hasher.update(fin.read(65000))
    return str(hasher.hexdigest())


#################################
#  File and Directory Handling  #
#################################


def get_test_dir():
    """Try to control for run location."""
    global ROOT_DIR
    global TEST_DIR
    if os.path.isfile('LICENSE') and os.path.isdir('ASEr'):
        ROOT_DIR = os.path.abspath('.')
    elif os.path.isfile('../LICENSE') and os.path.isdir('../ASEr'):
        ROOT_DIR = os.path.abspath('..')
    else:
        raise Exception('Run from project root')
    TEST_DIR = os.path.join(ROOT_DIR, 'tests', 'testdir')
    if not os.path.isdir(TEST_DIR):
        os.makedirs(TEST_DIR)


def get_file(infile):
    """Return path to a file, download if necessary."""
    path = os.path.join(TEST_DIR, infile)
    if os.path.isfile(path) and quick_hash(path) == FILES[infile]['short_hash']:
        return path
    else:
        return download_file(infile)
    raise Exception('Could not get file {}'.format(infile))


def download_file(infile):
    """Download the reference file."""
    path = os.path.join(TEST_DIR, infile)
    with open(path, 'wb') as fout:
        fout.write(urlopen(FILES[infile]['url']).read())
    if not quick_hash(path) == FILES[infile]['short_hash']:
        raise Exception('Downloaded file {} has incorrect hash'.format(path))
    return path


###############################################################################
#                               Test Functions                                #
###############################################################################


###################
#  Core Pipeline  #
###################


def test_countsnpase():
    """Run CountSNPASE.py with a random seed and compare to known hash."""
    # Make sure test dir is set and exists
    get_test_dir()

    # Build the command to run
    tmp_dir  = os.path.join(TEST_DIR, 'countsnpase_tmp')
    if not os.path.isdir(tmp_dir):
        os.makedirs(tmp_dir)
    command = """\
        python {root}/bin/CountSNPASE.py --mode single \
            --random-seed 0 \
            --reads {ase} \
            --snps {variants} \
            --prefix {tmpdir}/\
        """.format(root=ROOT_DIR, tmpdir=tmp_dir,
                   ase=get_file('ase.bam'), variants=get_file('variants.bed'))

    # Run the command and check success
    retcode, stdout, stderr = run.cmd(command)
    if not retcode == 0:
        sys.stderr.write('CODE: {}\nSTDOUT:\n{}\nSTDERR:\n{}\n'.format(
            retcode, stdout, stderr))
        raise Exception('CountSNPASE.py failed')

    # Check hash of the resulting file
    outfile  = os.path.join(tmp_dir, '_SNP_COUNTS.txt')
    snpfile  = os.path.join(TEST_DIR, 'SNP_COUNTS.txt')
    snp_hash = hash_file(outfile)
    assert snp_hash == FILES['SNP_COUNTS.txt']['hash']

    # Remove tmp files
    os.rename(outfile, snpfile)
    for f in os.listdir(tmp_dir):
        os.remove(os.path.join(tmp_dir, f))
    os.removedirs(tmp_dir)


def test_getgenease():
    """Run GetGeneASE.py and compare to known hash."""
    # Make sure test dir is set and exists
    get_test_dir()

    # Set the SNP file from the CountSNPASE step
    snpfile  = os.path.join(TEST_DIR, 'SNP_COUNTS.txt')
    outfile  = os.path.join(TEST_DIR, 'gene_ase.tsv')

    # Build the command to run
    tmp_dir  = os.path.join(TEST_DIR, 'countsnpase_tmp')
    if not os.path.isdir(tmp_dir):
        os.makedirs(tmp_dir)
    command = """\
        python bin/GetGeneASE.py \
            --snpcounts {snpfile} \
            --phasedsnps {variants} \
            --gff {gtf} \
            -o {tsvfile} \
            --writephasedsnps\
        """.format(snpfile=snpfile, variants=get_file('variants.bed'),
                   gtf=get_file('ref.gtf'), tsvfile=outfile)

    # Run the command and check success
    retcode, stdout, stderr = run.cmd(command)
    if not retcode == 0:
        sys.stderr.write('CODE: {}\nSTDOUT:\n{}\nSTDERR:\n{}\n'.format(
            retcode, stdout, stderr))
        raise Exception('GetGeneASE.py failed')

    # Check hash of the resulting file
    snp_hash = hash_file(outfile)
    assert snp_hash == FILES['gene_ase.tsv']['hash']

if __name__ == "__main__":
    test_countsnpase()
    test_getgenease()
    print("All tests successful!")
