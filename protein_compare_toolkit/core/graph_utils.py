'''
Modules for drawing and saving sequence logos.
'''
from Bio.Align import MultipleSeqAlignment
from pandas import DataFrame, RangeIndex

from matplotlib.figure import Figure
import matplotlib.pyplot as plt
import numpy as np
import logomaker as lm

from protein_compare_toolkit.core.stats_utils import info_content, select_diff_index, aln_dist

AA_SYMBOL = "ACDEFGHIKLMNPQRSTVWY"
AA_COUNT = 20

def plot_logo(
        aln1:MultipleSeqAlignment,
        aln2: MultipleSeqAlignment,
        start: int,
        end: int,
        plot_type: str
    ) -> Figure:
    '''A function that draws sequence logo.'''
    print(plot_type)
    if plot_type == "info":
        df1, df2 = info_logo_df(aln1, aln2)
        ylabel = "information content (bits)"
        ylim = (-4.33, 4.33)
    elif plot_type == "sdi":
        df1, df2 = sdi_logo_df(aln1, aln2)
        ylabel = "Selection-differentiation index"
        ylim = (-1., 1.)
    else:
        raise ValueError("Invalid Plot type.")

    # Reindex for plotting purpose.
    df1.index = RangeIndex(start, end)
    df2.index = RangeIndex(start, end)

    # Plotting logic.
    # Dynamically set figure width.
    width_per_residue = 0.5
    n = end - start + 1
    fig_width = max(4, n * width_per_residue)  # Set a minimum width for short windows.

    # Plot twice; one on the top, one on the bottom.
    fig, ax = plt.subplots(figsize=(fig_width, 4))
    lm.Logo(df=df1, ax=ax, color_scheme='skylign_protein')
    lm.Logo(df=df2, ax=ax, color_scheme="skylign_protein", fade_below=.5)

    # Adjust appearance.
    ax.axhline(0, color='black', linewidth=1)
    ax.set_xlim(start-.5, end+.5)
    ax.set_ylim(ylim)
    ax.set_xticks(np.arange(start, end+1, 1))
    ax.set_xlabel("Residue position")
    ax.set_ylabel(ylabel)

    fig.tight_layout()
    return fig


def info_logo_df(
        aln1: MultipleSeqAlignment,
        aln2: MultipleSeqAlignment) -> tuple[DataFrame, DataFrame]:
    '''A function that returns df needed to draw seq logo based on information content.'''
    info1 = info_content(aln1)
    info2 = info_content(aln2)
    dist1 = aln_dist(aln1)
    dist2 = aln_dist(aln2)

    df1 =  DataFrame(info1[:, np.newaxis] * dist1, columns=list(AA_SYMBOL))
    df2 = -DataFrame(info2[:, np.newaxis] * dist2, columns=list(AA_SYMBOL))
    return df1, df2


def sdi_logo_df(aln1: MultipleSeqAlignment, aln2: MultipleSeqAlignment):
    '''Make a pandas dataframe that the LogoMaker package can use to draw sequence logo with.'''
    sdi1, sdi2 = select_diff_index(aln1, aln2)
    dist1 = aln_dist(aln1)
    dist2 = aln_dist(aln2)

    # Total height represents sdi value;
    # Symbol height represents relative frequency.
    df1 =  DataFrame(sdi1[:, np.newaxis] * dist1, columns=list(AA_SYMBOL))
    df2 = -DataFrame(sdi2[:, np.newaxis] * dist2, columns=list(AA_SYMBOL))

    return df1, df2
