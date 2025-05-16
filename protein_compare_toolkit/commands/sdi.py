from pathlib import Path
from enum import Enum
from typing import Optional

import typer
from typing_extensions import Annotated
import pandas as pd

from protein_compare_toolkit.core.align_parser import read_alignment
from protein_compare_toolkit.core.graph_utils import plot_sdi_logo
from protein_compare_toolkit.core.stats_utils import select_diff_index, alignment_to_consensus

class SortBy(str, Enum):
    '''Enum for the --sort-by option of the rank command.'''
    aln1 = "aln1",
    aln2 = "aln2",
    avrg = "avrg"


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
    c_seq1 = alignment_to_consensus(sliced1)
    c_seq2 = alignment_to_consensus(sliced2)
    avg_sdi = (sdi1 + sdi2) / 2
    pos = list(range(start, end + 1))

    col = [("Position", ""),
           ("Alignment 1", "SDI"),
           ("Alignment 1", "Identity"),
           ("Alignment 2", "SDI"),
           ("Alignment 2", "Identity")]

    df = pd.DataFrame({
        "Position": pos,
        ("Alignment 1", "SDI"): sdi1,
        ("Alignment 1", "Identity"): c_seq1,
        ("Alignment 2", "SDI"): sdi2,
        ("Alignment 2", "Identity"): c_seq2
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
        typer.echo(df.to_string(index=False, justify="center"))
    else:
        df.to_csv(save_as, index=False)
        typer.echo(f"Ranking saved to {str(save_as)}")
