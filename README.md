# Protein Compare Toolkit CLI

Protein Compare Toolkit is a command-line tool to analyze and visualize the **Selection-Differentiation Index (SDI)** between two protein sequence alignments.

## Installation

```bash
git clone https://github.com/PhantomRabbit/protein_compare_toolkit
```

Head to project directory.

```bash
pip install -e .
```

## Usage

### SDI RANK

To compare two protein alignments using the selection-differentiation index and generate a ranked table, use the following command.

```bash
protein-compare-toolkit sdi rank <alignment1>.aln <alignment2>.aln [<start>] [<end>] [<top>] [--by-first|--by-second] [--save-as <rank>.csv]
```

#### Arguments

- `alignment1`. Name of the first multiple sequence alignment file.
- `alignment2`. Name of the second multiple sequence alignment file.
- `start`. The start position on the alignment to be graphed (1-index and inclusive). Default to the start of the alignment.
- `end`. The end position. Default to the end of the alignment.
- `top`. How many high-score positions should the ranking return. Default is 10.

#### Options

- `--by-first|--by-second`. The alignment from which the SDI value should be ranked. Default to the first alignment.
- `--peek|--save`. Whether the ranking should be displayed in the terminal, or be saved as a CSV file. Default to `--peek`.
- `--save-as`. The file that the output table should be saved as when the output mode is `--save`.

### SDI LOGO

To compare two protein sequence alignments and generate a sequence logo using the selection-differentiation index, use the following command:

```bash
protein-compare-toolkit sdi logo <alignment1>.aln <alignment2>.aln [<start>] [<end>] [--save-as <logo>.png]
```

#### Arguments

- `alignment1`. Name of the first multiple sequence alignment file.
- `alignment2`. Name of the second multiple sequence alignment file.
- `start`. The start position. Default to the start of the alignment.
- `end`. The end position. Default to the end of the alignment.

#### Options

- `--save-as`. The file that the output graph should be saved as.

## License

