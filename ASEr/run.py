"""
File management and execution functions.

============================================================================

        AUTHOR: Michael D Dacre, mike.dacre@gmail.com
  ORGANIZATION: Stanford University
       LICENSE: MIT License, property of Stanford, use as you wish
       VERSION: 0.1
       CREATED: 2016-02-11 16:03
 Last modified: 2016-03-14 12:42

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
    changes."""
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


def which(program):
    """Replicate's the UNIX which command.

    Taken verbatim from:
        stackoverflow.com/questions/377017/test-if-executable-exists-in-python

    :program: Name of executable to test.
    :returns: Path to the program or None on failure.
    """
    def is_exe(fpath):
        """Return True is fpath is executable."""
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

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
