'''Commands that plots a sequence logo using different metrics.'''
from pathlib import Path
from enum import Enum

import typer
from typing_extensions import Annotated

from protein_compare_toolkit.core.align_parser import read_alignment, slice_alignment
from protein_compare_toolkit.core import graph_utils

class Metric(str, Enum):
    '''Metric used for plotting the sequence logo.'''
    INFO = "info"
    SDI = "sdi"

app = typer.Typer()

@app.command()
def logo(
    metric: Annotated[
        Metric,
        typer.Argument(help="Choose the metric you want to use for plotting. ('info' or 'sdi')")
    ],
    file1: Annotated[
        Path,
        typer.Argument(exists=True, file_okay=True, readable=True)],
    file2: Annotated[
        Path,
        typer.Argument(exists=True, file_okay=True, readable=True)],
    start: Annotated[
        int,
        typer.Argument(
            min=1,
            clamp=True,
            show_default=False,
            help="Start position of the graphing range (inclusive)"
        )
    ] = -1,
    end: Annotated[
        int,
        typer.Argument(
            min=1,
            clamp=True,
            show_default=False,
            help="End position of the graphing range (inclusive)"
        )
    ] = -1,
    output: Annotated[
        Path,
        typer.Argument(
            file_okay=True,
            writable=True,
            help="Output image file"
        )
    ] = Path("logo.png")
):
    '''
    The logo let the user draw a sequence logo based on selected metrics.

    Usage: protein-compare-toolkit logo ["info"/"sdi"] [file1] [file2] [start] [end] [output]
    '''

    aln1 = read_alignment(file1)
    aln2 = read_alignment(file2)

    aln_len = len(aln1[0])
    if aln_len != len(aln2[0]):
        raise ValueError(
            f"Alignment lengths differ: {len(aln1[0])} vs {len(aln2[0])}"
        )
    # Handle out of bound and impossible indices for the range.
    if start == -1 or start > aln_len:
        start = 1
    if end == -1 or end > aln_len:
        start = aln_len
    if start > end:
        start = 1
        end = aln_len

    sliced1 = slice_alignment(aln1, start, end)
    sliced2 = slice_alignment(aln2, start, end)

    fig = graph_utils.plot_logo(sliced1, sliced2, start, end, metric)
    fig.savefig(output)
    typer.echo(f"Graph saved to {output}.")
