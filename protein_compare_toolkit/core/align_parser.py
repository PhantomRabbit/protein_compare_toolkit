'''
This module contains I/O stuff.

Not much else to say (yet).
'''
from pathlib import Path

from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment

def read_alignment(aln_path: Path) -> MultipleSeqAlignment:
    '''Reads in a clustal multiple sequence alignment file.'''
    aln = AlignIO.read(aln_path, "clustal")
    return aln

def slice_alignment(aln: MultipleSeqAlignment, start: int, end: int) -> MultipleSeqAlignment:
    '''Slice a multiple sequence alignment object using 1-base inclusive indices.'''
    sliced = aln[:, start:end+1]
    assert isinstance(sliced, MultipleSeqAlignment)
    return sliced
