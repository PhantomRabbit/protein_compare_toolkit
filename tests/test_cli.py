from pathlib import Path

from typer.testing import CliRunner
from protein_compare_toolkit.commands.sdi import app

runner = CliRunner()

def test_cli_generates_output(tmp_path):
    # Locate test alignment files relative to this test file
    test_dir = Path(__file__).parent
    data_dir = test_dir / "data"

    file1 = data_dir / "fesii_align.aln"
    file2 = data_dir / "non_fesii_align.aln"
    output_file = tmp_path / "cli_output.png"

    result = runner.invoke(
        app,
        [
            # "sequences",
            str(file1),
            str(file2),
            "--start", "75",
            "--end", "125",
            "--output", str(output_file)
        ]
    )

    assert result.exit_code == 0, result.stdout
    assert output_file.exists()
    assert output_file.stat().st_size > 0
