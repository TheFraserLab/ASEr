"""
Get ASE counts from BAMs or raw fastq data.

============================================================================

        AUTHOR: Michael D Dacre, mike.dacre@gmail.com
       CREATOR: Carlo Artieri
       CREATOR: Tomas Babak
  ORGANIZATION: Stanford University
       LICENSE: MIT License, property of Stanford, use as you wish
       VERSION: 0.3.0
       CREATED: 2016-25-15 17:03
 Last modified: 2016-03-24 16:19

============================================================================
"""
from . import snps
from . import plink

__all__ = ['snps', 'plink']
