'''
Statistics functions for protein sequence comparison.

This module includes the main function that calculates the detailed matrix to draw sequence logo.
It also includes functions that calculates the two key metrics, Jensen-Shannon Distance and 
information content, and associated utility functions.
'''
from Bio.Align import MultipleSeqAlignment
from scipy import spatial, stats
import numpy as np

AA_SYMBOL = "ACDEFGHIKLMNPQRSTVWY"
AA_COUNT = 20

def select_diff_index(aln1: MultipleSeqAlignment, aln2: MultipleSeqAlignment):
    '''
    Calculate the selection-differentiation index.

    The selection-differentiation index is a metric that quantifies the divergence of two related
    families of proteins on a residue sequence basis. The formula for the index at the nth position
    from family A is:
    
    SDI_A_n = I_A_n * JSD_n / log2(20)
    
    where I_A_n is the information content at the nth position in family A, a metric of selective
    pressure; and JSD_n is the statistical distance between the distributions of residues in the two 
    alignments at this position, which indicates differentiation between the two families. Then the
    metric is normalised to 1, dividing by the maximum possible information content of any residue
    distribution.

    This composite metric highlights residues that are under differentiating selective pressure, 
    which are likely to be responsible for functional difference.
    '''
    info1 = info_content(aln1)
    info2 = info_content(aln2)
    jsd = jensen_shannon_distance(aln1, aln2)

    sdi1 = info1 * jsd / np.log2(AA_COUNT) #Normalisation.
    sdi2 = info2 * jsd / np.log2(AA_COUNT)

    return sdi1, sdi2


def info_content(alignment: MultipleSeqAlignment):
    '''
    Calculates the information content of positions in a alignment.
    
    This the information content of a residue distribution in a protein sequence when compared to
    background distribution (uniform). It is a metric for measuring the conservation of a particular
    residue, hinting at functional significance due to selective pressure.
    '''
    if len(alignment) == 0:
        raise ValueError("Alignment must contain at least one sequence.")

    seq_count = len(alignment)
    alignment_len = len(alignment[0])
    r0 = np.log2(20)
    r = np.zeros(alignment_len)
    f = alignment_to_distribution(alignment)
    e = 1 / np.log(2) * (20 - 1) / (2 * seq_count) #Error correction.

    for i in range(alignment_len): #Formula reference: https://en.wikipedia.org/wiki/Sequence_logo
        h = stats.entropy(f[i])
        r[i] = max(0, r0 - (h + e))  # Clamp negatives to 0

    return r


def jensen_shannon_distance(aln1: MultipleSeqAlignment, aln2: MultipleSeqAlignment):
    '''
    The statistical distance between two distributions of residues in two matching alignments.
    
    This metric is selected for the explicit reason that it is symmetric, i.e. the distance from A
    to B is the same as that from B to A (unlike statistical divergence). Another reason is that
    the metric is bounded between 0 and 1, which is easier to handle and understand.
    '''
    validate_alignment_lengths(aln1, aln2)

    aln_len = len(aln1[0])
    jsd = np.zeros(aln_len)

    p = alignment_to_distribution(aln1)
    q = alignment_to_distribution(aln2)

    for k in range(aln_len):
        jsd[k] = spatial.distance.jensenshannon(p[k], q[k], base=2)

    return jsd


def alignment_to_distribution(aln: MultipleSeqAlignment):
    '''Normalise counts of residues in a alignment so the sum is 1 for all positions.'''
    aln_len = len(aln[0])
    dist = np.zeros((aln_len, AA_COUNT))

    for i in range(aln_len):
        d = [aln[:, i].count(aa) for aa in AA_SYMBOL]
        s = np.sum(d)
        if s == 0:
            d = np.full(AA_COUNT, 1 / AA_COUNT)
            # return a normal dist for gaps for better error handling.
        else:
            d = np.array(d) / s
        dist[i] = d

    return dist


def validate_alignment_lengths(aln1: MultipleSeqAlignment, aln2: MultipleSeqAlignment):
    '''Raises ValueError if alignments do not have the same length.'''
    if len(aln1[0]) != len(aln2[0]):
        raise ValueError(
            f"Alignment lengths differ: {len(aln1[0])} vs {len(aln2[0])}"
        )
