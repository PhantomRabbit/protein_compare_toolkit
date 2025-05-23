from pathlib import Path

import pandas as pd
import typer
from typing_extensions import Annotated

from protein_compare_toolkit.core.align_parser import read_alignment, slice_alignment
from protein_compare_toolkit.core.stats_utils import consensus_seq, select_diff_index

app = typer.Typer()

@app.command()
def rank(
    file1: Annotated[Path, typer.Argument(exists=True, file_okay=True, readable=True)],
    file2: Annotated[Path, typer.Argument(exists=True, file_okay=True, readable=True)],
    start: Annotated[int, typer.Argument(
        help="Start position of the graphing range (inclusive)",
        show_default=False
        )] = -1,
    end: Annotated[int, typer.Argument(
        help="End position of the graphing range (inclusive)",
        show_default=False
        )] = -1,
    top: Annotated[int, typer.Argument(
        help="Return only the top elements from the ranking",
        show_default=False
    )] = 25,
    sort_by: Annotated[bool, typer.Option("--by-1st/--by-2nd")] = True,
    output_mode: Annotated[bool, typer.Option("--peek/--save")] = True,
    save_as: Annotated[Path, typer.Option(
        file_okay=True,
        writable=True,
        help="Path to save CSV file when the output mode is 'save'."
        )] = Path("sdi_rank.csv")
):
    '''
    The rank commands ranks alignment positions based on their SDI value.

    Usage: protein-compare-toolkit sdi rank <file1> <file2> <start> <end> 
    --top <n> [--peek/--save] --save-as <path>.
    '''

    aln1 = read_alignment(file1)
    aln2 = read_alignment(file2)
    aln_len = len(aln1[0]) # Aligned sequences should all have equal lengths.

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

    sdi1, sdi2 = select_diff_index(sliced1, sliced2)
    c_seq1 = consensus_seq(sliced1)
    c_seq2 = consensus_seq(sliced2)

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
        "Position": list(range(start, end)),
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
        by = ("Alignment 1", "SDI")
    else:
        by = ("Alignment 2", "SDI")

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
