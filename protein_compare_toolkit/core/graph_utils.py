'''
Modules for drawing and saving sequence logos.
'''
from Bio.Align import MultipleSeqAlignment
from pandas import DataFrame, RangeIndex

import matplotlib.pyplot as plt
import numpy as np
import logomaker as lm

from protein_compare_toolkit.core.stats_utils import alignment_to_distribution, select_diff_index

AA_SYMBOL = "ACDEFGHIKLMNPQRSTVWY"
AA_COUNT = 20

def plot_sdi_logo(aln1: MultipleSeqAlignment, aln2: MultipleSeqAlignment,
                  start: int, end: int, output_file: str):
    '''
    Draw a sequence logo using the selection-differentiation index of the two families of proteins.

    Takes a pair of sequence, a start and an end position, and a target output file.
    '''
    df1, df2 = sdi_logo_matrix(aln1, aln2, start, end)
    df2 = -df2 # Flip bottom logo values to be drawn on the bottom.
    # Since Logomaker doesn't support using repeated glyph in the same logo, this code forces it to
    # draw two logos with positive and negative height respectively on the same axes to achieve the
    # effct.

    # Dynamically set figure width
    width_per_residue = 0.5  # Adjust this to control column width
    n_residues = end - start + 1
    fig_width = max(4, n_residues * width_per_residue)  # Set a minimum width for short windows

    fig, ax = plt.subplots(figsize=(fig_width, 4))
    lm.Logo(df=df1, ax=ax, color_scheme='skylign_protein')
    lm.Logo(df=df2, ax=ax, color_scheme="skylign_protein", fade_below=.5)
    # Adjust appearance.
    ax.axhline(0, color='black', linewidth=1)
    ax.set_ylim(-1, 1)
    ax.set_xlim(start - 1, end)
    ax.set_xlabel("Residue position")
    ax.set_ylabel("Selection-differentiation index")

    fig.tight_layout()
    fig.savefig(output_file)
    plt.close(fig)


def sdi_logo_matrix(aln1: MultipleSeqAlignment, aln2: MultipleSeqAlignment, start: int, end: int):
    '''Make a pandas dataframe that the LogoMaker package can use to draw sequence logo with.'''
    sdi1, sdi2 = select_diff_index(aln1, aln2)
    dist1 = alignment_to_distribution(aln1)
    dist2 = alignment_to_distribution(aln2)

    df1 = DataFrame(sdi1[:, np.newaxis] * dist1, RangeIndex(start - 1, end), list(AA_SYMBOL))
    df2 = DataFrame(sdi2[:, np.newaxis] * dist2, RangeIndex(start - 1, end), list(AA_SYMBOL))
    # Total height represents sdi value;
    # Symbol height represents relative frequency.

    return df1, df2
