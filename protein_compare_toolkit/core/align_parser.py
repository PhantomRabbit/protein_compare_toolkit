'''
This module contains I/O stuff.

Not much else to say (yet).
'''
from Bio import Align, AlignIO

def read_alignment(file_path: str) -> Align.MultipleSeqAlignment:
    '''Reads in a clustal multiple sequence alignment file.'''
    alignment = AlignIO.read(file_path, "clustal")
    return alignment
