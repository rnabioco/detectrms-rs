
# Detect RMS (RNA Modification Sites)

This tool processes BAM alignment files from Nanopore long read sequencing and using the references fasta file, computes the mismatch frequency of all bases among other stats.

## Features
- Reading and processing BAM files using `rust_htslib`.
- Extracting sequences and their information from FASTA files using `bio::io::fasta`.
- Calculating statistical measures such as mean, median, and standard deviation of quality scores.
- Generating pileup information and computing mismatch, insertion, and deletion ratios.
- Writing the processed data to a TSV file for easy analysis and visualization.

## Installation

To install `detectrms`, you need to have Rust and Cargo installed on your system. If you don't have Rust installed, follow the instructions on [Rust's official website](https://www.rust-lang.org/tools/install).

```bash
git clone https://github.com/rnabioco/detectrms-rs

cd detectrms-rs

cargo install --path .
```

## Usage
After building the project, you can run the tool using the following command:

```
detectrms --bam <path-to-bam-file> --fasta <path-to-fasta-file> --output <path-to-output-tsv>
```

### Command Line Arguments
- `--bam`: Path to the input BAM file.
- `--fasta`: Path to the input FASTA file.
- `--output`: Path where the output TSV file will be saved.

## Output Format
The output is a TSV file with the following columns:
- Reference ID
- Position
- Base
- Strand
- Coverage
- Mean Quality Score
- Median Quality Score
- Standard Deviation of Quality Scores
- Mismatch Ratio
- Insertion Ratio
- Deletion Ratio
- Counts of A, C, G, T bases

## Note
Idea and formmated output is derived from nanoRMS. <https://github.com/novoalab/nanoRMS>
