# amplipore

Amplicon analysis of Oxford Nanopore 16S V1V9 sequences.

## Requirements

- Docker Engine

## Installation

Clone the repository or download the source code.

Build image with `make build`.

## Usage

All commands need to be run inside the container. Use `make bash` for this.

The directories `app` (sourcecode) and `data` (put files here) are volume-mounted into the container.

## Building database

Use `database.py` to download the database, which is a reformatted version of the [EMU](https://github.com/treangenlab/emu) database.

`python database.py --outdir <outdir>`

## Classification

The reads are preprocessed with [fastq_rs](https://github.com/OscarAspelin95/fastq_rs), clustered with [USEARCH12](https://github.com/rcedgar/usearch12) and classified with [sintax_rs](https://github.com/OscarAspelin95/sintax_rs).

`python main.py --fastq <reads.fastq.gz> [...] --database <db.fasta> --outdir <outdir>`

Optional arguments:

<pre>
<b>-s/--sintax_threshold</b> [0.80] - Threshold for assigning a taxonomic level.
</pre>

## Output Files

All result files are generated in the `outdir` directory.
