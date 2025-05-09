# Protein Compare Toolkit CLI

Protein Compare is a command-line tool to analyze and visualize the **Selection-Differentiation Index (SDI)** between two protein sequence alignments.

## Installation

## Usage

To compare two protein sequence alignments and generate a sequence logo, use the following command:

```bash
protein-compare-toolkit sdi logo <file1> <file2> --start <start_position> --end <end_position> --output <output_file>
```

### Options:

- <file1>: Path to the first protein alignment (e.g., alignment1.aln).

- <file2>: Path to the second protein alignment (e.g., alignment2.aln).

- --start: The start position for the alignment range (inclusive).

- --end: The end position for the alignment range (inclusive).

- --output: Path to save the resulting image (default: comparison.png).

## License

