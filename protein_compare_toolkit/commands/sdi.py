import typer
from protein_compare_toolkit.core.align_parser import read_alignment
from protein_compare_toolkit.core.graph_utils import plot_sdi_logo

app = typer.Typer()

@app.command()
def logo(
    file1: str,
    file2: str,
    start: int = typer.Option(..., help="Start position of the range"),
    end: int = typer.Option(..., help="End position of the range"),
    output: str = typer.Option("sdi_logo.png", help="Output image file")
):
    aln1 = read_alignment(file1)
    aln2 = read_alignment(file2)

    sliced1 = aln1[:, start-1:end]
    sliced2 = aln2[:, start-1:end]

    plot_sdi_logo(sliced1, sliced2, start, end, output)
    typer.echo(f"Graph saved to {output}.")
