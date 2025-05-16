'''
This test contains validation for the function of graphing utilities.
'''
from pathlib import Path
from protein_compare_toolkit.core.align_parser import read_alignment
from protein_compare_toolkit.core.graph_utils import plot_sdi_logo

def test_plot_sdi_logo_creates_file(tmp_path):
    '''Test if the plotting function correctly produces a file.'''
    # Get the directory where THIS test file lives
    test_dir = Path(__file__).parent
    data_dir = test_dir / "data"

    aln1 = read_alignment(data_dir / "fesii_align.aln")
    aln2 = read_alignment(data_dir / "non_fesii_align.aln")

    aln1 = aln1[:, 75 - 1:125]
    aln2 = aln2[:, 75 - 1:125]

    output_file = tmp_path / "sdi_logo_debug.png"
    plot_sdi_logo(aln1, aln2, 75, 125, str(output_file))

    assert output_file.exists()
    assert output_file.stat().st_size > 0

