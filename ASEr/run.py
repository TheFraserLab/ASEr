"""
File management and execution functions.

============================================================================

        AUTHOR: Michael D Dacre, mike.dacre@gmail.com
  ORGANIZATION: Stanford University
       LICENSE: MIT License, property of Stanford, use as you wish
       VERSION: 0.1
       CREATED: 2016-02-11 16:03
 Last modified: 2016-03-21 13:54

   DESCRIPTION: Run commands with run_cmd, search the PATH with which.

         USAGE: Import as a module or run as a script

============================================================================
"""
import os
import gzip
import bz2
from subprocess import Popen
from subprocess import PIPE

from . import logme

__all__ = ['run_cmd', 'which', 'open_zipped']


def open_zipped(infile, mode='r'):
    """Open a regular, gzipped, or bz2 file.

    Returns text mode file handle.

    If infile is a file handle or text device, it is returned without
    changes.
    """
    mode   = mode[0] + 't'
    if hasattr(infile, 'write'):
        return infile
    if isinstance(infile, str):
        if infile.endswith('.gz'):
            return gzip.open(infile, mode)
        if infile.endswith('.bz2'):
            if hasattr(bz2, 'open'):
                return bz2.open(infile, mode)
            else:
                return bz2.BZ2File(infile, mode)
        return open(infile, mode)


def run_cmd(cmd, args):
    """Run command and return status, output, stderr.

    cmd:  Path to executable.
    args: Tuple of arguments.
    """
    args = (cmd,) + args
    pp = Popen(args, shell=False, universal_newlines=True,
               stdout=PIPE, stderr=PIPE)
    out, err = pp.communicate()
    code = pp.returncode
    if out[-1:] == '\n':
        out = out[:-1]
    if err[-1:] == '\n':
        err = err[:-1]
    return code, out, err


def is_exe(fpath):
    """Return True is fpath is executable."""
    return os.path.isfile(fpath) and os.access(fpath, os.X_OK)


def which(program):
    """Replicate the UNIX which command.

    Taken verbatim from:
        stackoverflow.com/questions/377017/test-if-executable-exists-in-python

    :program: Name of executable to test.
    :returns: Path to the program or None on failure.
    """
    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return os.path.abspath(program)
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return os.path.abspath(exe_file)

    return None


def file_type(infile):
    """Return file type after stripping gz or bz2."""
    name_parts = infile.split('.')
    if name_parts[-1] == 'gz' or name_parts[-1] == 'bz2':
        name_parts.pop()
    return name_parts[-1]


def is_file_type(infile, types):
    """Return True if infile is one of types.

    :infile:  Any file name
    :types:   String or list/tuple of strings (e.g ['bed', 'gtf'])
    :returns: True or False

    """
    if hasattr(infile, 'write'):
        return False
    if isinstance(types, str):
        types = [types]
    if not isinstance(types, (list, tuple)):
        raise Exception('types must be string list or tuple')
    for typ in types:
        if file_type(infile) == typ:
            return True
    return False


def write_iterable(iterable, outfile):
    """Write all elements of iterable to outfile."""
    with open_zipped(outfile, 'w') as fout:
        fout.write('\n'.join(iterable))


def split_file(infile, parts, outpath='', keep_header=True):
    """Split a file in parts parts and return a list of paths.

    NOTE: Linux specific (uses wc).

    :outpath:     The directory to save the split files.
    :keep_header: Add the header line to the top of every file.

    """
    # Determine how many reads will be in each split sam file.
    logme.log('Getting line count', 'debug')
    num_lines = int(os.popen(
        'wc -l ' + infile + ' | awk \'{print $1}\'').read())
    num_lines   = int(int(num_lines)/int(parts)) + 1

    # Subset the file into X number of jobs, maintain extension
    cnt      = 0
    currjob  = 1
    suffix   = '.split_' + str(currjob).zfill(4) + '.' + infile.split('.')[-1]
    file_path, file_name = os.path.split(infile)
    run_file = os.path.join(outpath, file_name + suffix)
    outfiles = [run_file]

    # Actually split the file
    logme.log('Splitting file', 'debug')
    with open(infile) as fin:
        header = fin.readline() if keep_header else ''
        split_file = open(run_file, 'w')
        split_file.write(header)
        for line in fin:
            cnt += 1
            if cnt < num_lines:
                split_file.write(line)
            elif cnt == num_lines:
                split_file.write(line)
                split_file.close()
                currjob += 1
                suffix = '.split_' + str(currjob).zfill(4) + '.' + \
                    infile.split('.')[-1]
                run_file = os.path.join(outpath, file_name + suffix)
                split_file = open(run_file, 'w')
                outfiles.append(run_file)
                split_file.write(header)
                cnt = 0
        split_file.close()
    return tuple(outfiles)

