from pathlib import Path

from typer.testing import CliRunner
from protein_compare_toolkit.commands.sdi import app

runner = CliRunner()

def test_logo_generates_output(tmp_path):
    '''Test if the logo command correctly generate an output file.'''
    # Locate test alignment files relative to this test file
    test_dir = Path(__file__).parent
    data_dir = test_dir / "data"

    file1 = data_dir / "fesii_align.aln"
    file2 = data_dir / "non_fesii_align.aln"
    output_file = tmp_path / "cli_output.png"

    result = runner.invoke(
        app,
        [
            "logo",
            str(file1),
            str(file2),
            "75",
            "125",
            "--save-as", str(output_file)
        ]
    )

    assert result.exit_code == 0, result.stdout
    assert output_file.exists()
    assert output_file.stat().st_size > 0


def test_rank_generates_output(tmp_path):
    '''Test if the rank command correctly generate an output file.'''
    test_dir = Path(__file__).parent
    data_dir = test_dir / "data"

    file1 = data_dir / "fesii_align.aln"
    file2 = data_dir / "non_fesii_align.aln"
    output_file = tmp_path / "cli_output.csv"

    result = runner.invoke(
        app,
        [
            "rank",
            str(file1),
            str(file2),
            "75",
            "125",
            "--save",
            "--save-as", str(output_file)
        ]
    )

    assert result.exit_code == 0, result.stdout
    assert output_file.exists()
    assert output_file.stat().st_size > 0

