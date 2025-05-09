'''
Contains test files for testing the functionality of the project.
'''
from pathlib import Path
from protein_compare_toolkit.core.align_parser import read_alignment


def test_read_alignment_returns_alignment_object():
    '''test if the read_alignment function correctly reads alignment from test file.'''
    file_path = Path(__file__).parent / "data" / "fesii_align.aln"
    alignment = read_alignment(file_path)
    print(alignment)
    assert len(alignment) > 1
    assert alignment.get_alignment_length() > 0
