"""
Simple wrapper for some plink commands relevant to this pipleine.

============================================================================

        AUTHOR: Michael D Dacre, mike.dacre@gmail.com
  ORGANIZATION: Stanford University
       LICENSE: MIT License, property of Stanford, use as you wish
       VERSION: 0.1
       CREATED: 2016-31-11 14:03
 Last modified: 2016-03-18 16:11

   DESCRIPTION:

         USAGE: Import as a module or run as a script

============================================================================
"""
import os
import sys

# Us
from .run import run_cmd
from .run import which
from .run import open_zipped

# Logging
from . import logme
logme.MIN_LEVEL = 'info'

__all__ = ['plink', 'is_recodeAD', 'recodeAD']

###############################################################################
#                             Run Plink Commands                              #
###############################################################################


def plink(args, plink_exec=None, logfile=None):
    """Run plink with *args

    :args:       List, tuple, or string of arguments to pass to plink.
    :plink_exec: Location of plink, if not provided, PATH searched.
    :logfile:    A file to write log messages too, not required.
    :returns:    stdout, stderr

    """
    # Set logfile
    if logfile:
        logme.LOGFILE = logfile
    else:
        logfile = sys.stderr

    # Make sure plink exists
    plink_exec = os.path.abspath(plink) if plink_exec else which('plink')
    if not plink:
        raise PlinkError('Could not find plink in the PATH')
    if not os.path.isfile(plink_exec) or not os.access(plink_exec, os.X_OK):
        raise PlinkError('{} is not executable'.format(plink_exec))

    if isinstance(args, list):
        args = tuple(args)
    elif isinstance(args, str):
        args = tuple(args.split(' '))
    if not isinstance(args, tuple):
        raise PlinkError('args must be string, tuple, or list')

    # Actually run the command
    code, stdout, stderr = run_cmd(plink_exec, args)

    # Check the results
    if code is not 0:
        logme.log('Plink run failed with args: {}\n\n'.format(args),
                  level='error')
        sys.stderr.write('STDOUT:\n{}\n\nSTDERR:{}\n'.format(stdout, stderr))
        raise PlinkError('Plink failed with exit code {}'.format(code))

    return stdout, stderr


###############################################################################
#                             Tests of file type                              #
###############################################################################


def is_recodeAD(infile):
    """Check a .raw file. Return True if is recodeAD file, False otherwise."""
    with open_zipped(infile) as fin:
        return True if fin.read(500).split(' ')[7].endswith('HET') else False


###############################################################################
#                            Conversion Functions                             #
###############################################################################


def recodeAD(infile, plink_exec=None):
    """Check if infile is in AD format and if not, recodeAD.

    :infile:     A plink prefix.
    :plink_exec: Location of plink, if not provided, PATH searched.
    :returns:    The path to the raw file.

    """
    infile    = get_root_name(infile)
    outfile   = infile + '.recodeAD'
    file_flag = get_file_flag(infile)
    if not file_flag:
        raise PlinkError("Coundn't find either a bed or ped file")
    try:
        plink((file_flag, infile, '--out', outfile, '--recodeAD'),
               plink_exec=plink_exec)
    except PlinkError:
        logme.log('Failed during recodeAD step', level='error')
        raise

    return os.path.abspath(outfile + '.raw')


###############################################################################
#                           Housekeeping Functions                            #
###############################################################################


def get_file_flag(infile):
    """Return bfile if bed exists of file if ped exists."""
    infile = get_root_name(infile)
    if os.path.isfile(infile + '.bed'):
        return '--bfile'
    elif os.path.isfile(infile + '.ped'):
        return '--file'
    else:
        return None


def get_root_name(infile):
    """Return root plink prefix for file with a known filetype."""
    file_endings = ['bed', 'bim', 'fam', 'map', 'ped', 'raw', 'tfam', 'tped']
    for ending in file_endings:
        if infile.endswith(ending):
            return '.'.join(infile.split('.')[:-1])
    return infile



class PlinkError(logme.LoggingException):

    """A Logging Excption Class."""

    pass
