'''
All the commands relevant to the use of Selection-Differentiation Index.

logo: plot a sequence logo within given range.
rank: rank alignment positions by their SDI value.
'''
from pathlib import Path
from typing import Optional

import typer
from typing_extensions import Annotated
import pandas as pd

from protein_compare_toolkit.core.align_parser import read_alignment
from protein_compare_toolkit.core.graph_utils import plot_sdi_logo
from protein_compare_toolkit.core.stats_utils import select_diff_index, consensus_seq

app = typer.Typer()

@app.command()
def logo(
    file1: Annotated[Path, typer.Argument(exists=True, file_okay=True, readable=True)],
    file2: Annotated[Path, typer.Argument(exists=True, file_okay=True, readable=True)],
    start: Annotated[
        int,
        typer.Argument(
            min=1,
            clamp=True,
            help="Start position of the graphing range (inclusive)"
            )
        ],
    end: Annotated[
        int,
        typer.Argument(
            min=1,
            clamp=True,
            help="End position of the graphing range (inclusive)"
            )
        ],
    save_as: Annotated[
        Optional[Path],
        typer.Option(
            file_okay=True,
            writable=True,
            help="Output image file"
            )
        ] = Path("logo.png"),
):
    '''
    The logo commands let the user draw a sequence logo based on selection-differentiation index.

    Usage: protein-compare-toolkit sdi logo <file1> <file2> <start> <end> --save-as <path>
    '''
    aln1 = read_alignment(file1)
    aln2 = read_alignment(file2)

    sliced1 = aln1[:, start-1:end]
    sliced2 = aln2[:, start-1:end]

    plot_sdi_logo(sliced1, sliced2, start, end, save_as)
    typer.echo(f"Graph saved to {str(save_as)}.")

@app.command()
def rank(
    file1: Annotated[Path, typer.Argument(exists=True, file_okay=True, readable=True)],
    file2: Annotated[Path, typer.Argument(exists=True, file_okay=True, readable=True)],
    start: Annotated[
        int,
        typer.Argument(
            help="Start position of the graphing range (inclusive)",
            show_default=False
            )
        ] = 1,
    end: Annotated[
        int,
        typer.Argument(
            help="End position of the graphing range (inclusive)",
            show_default=False
            )
        ] = -1,
    top: Annotated[
        int,
        typer.Option(
            help="Return only the top elements from the ranking",
            show_default=False
            )
        ] = 10,
    sort_by: Annotated[Optional[bool], typer.Option("--by-1st/--by-2nd")] = True,
    output_mode: Annotated[Optional[bool], typer.Option("--peek/--save")] = True,
    save_as: Annotated[
        Optional[Path],
        typer.Option(
            file_okay=True,
            writable=True,
            help="Path to save CSV file when the output mode is 'save'."
            )
        ] = Path("sdi_rank.csv")
):
    '''
    The rank commands ranks alignment positions based on their SDI value.

    Usage: protein-compare-toolkit sdi rank <file1> <file2> <start> <end> 
    --top <n> [--peek/--save] --save-as <path>.
    '''

    aln1 = read_alignment(file1)
    aln2 = read_alignment(file2)
    aln_len = len(aln1[0]) # Aligned sequences should all have equal lengths.

    # Check for invalid values.
    if start < 1 or start > aln_len:
        start = 1

    if end < 1 or end > aln_len:
        end = aln_len

    if end <= start:
        start = 1
        end = aln_len

    sliced1 = aln1[:, start-1:end]
    sliced2 = aln2[:, start-1:end]

    sdi1, sdi2 = select_diff_index(sliced1, sliced2)
    c_seq1 = consensus_seq(sliced1)
    c_seq2 = consensus_seq(sliced2)
    pos = list(range(start, end + 1))

    col = [("Position", ""),
           ("Alignment 1", "SDI"),
           ("Alignment 1", "Identity"),
           ("Alignment 1", "P"),
           ("Alignment 1", "Lower"),
           ("Alignment 1", "Upper"),
           ("Alignment 2", "SDI"),
           ("Alignment 2", "Identity"),
           ("Alignment 2", "P"),
           ("Alignment 2", "Lower"),
           ("Alignment 2", "Upper")
           ]

    df = pd.DataFrame({
        "Position": pos,
        ("Alignment 1", "SDI"): sdi1,
        ("Alignment 1", "Identity"): list(c_seq1["id"]),
        ("Alignment 1", "P"): list(c_seq1["p"]),
        ("Alignment 1", "Lower"): list(c_seq1["lower"]),
        ("Alignment 1", "Upper"): list(c_seq1["upper"]),
        ("Alignment 2", "SDI"): sdi2,
        ("Alignment 2", "Identity"): list(c_seq2["id"]),
        ("Alignment 2", "P"): list(c_seq2["p"]),
        ("Alignment 2", "Lower"): list(c_seq2["lower"]),
        ("Alignment 2", "Upper"): list(c_seq2["upper"])
    })

    if sort_by:
        by = [("Alignment 1", "SDI")]
    else:
        by = [("Alignment 2", "SDI")]

    df.columns = pd.MultiIndex.from_tuples(col)
    df = df.sort_values(by=by, ascending=False)

    if top > aln_len:
        typer.echo("Input 'top' value is invalid. Returning the entire ranking.")
        top = aln_len

    df = df[:top]

    if output_mode: # --peek
        typer.echo(df.to_string(index=False, justify="center", float_format='{:,.3f}'.format))
    else:
        df.to_csv(save_as, index=False)
        typer.echo(f"Ranking saved to {str(save_as)}")
