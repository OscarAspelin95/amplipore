# amplipore

16S V1-V9 amplicon classification pipeline for Oxford Nanopore reads.

Reads are preprocessed with [fastq_rs](https://github.com/OscarAspelin95/fastq_rs), clustered with [USEARCH12](https://github.com/rcedgar/usearch12), and classified with [sintax_rs](https://github.com/OscarAspelin95/sintax_rs). Optionally, unclassified ASVs are re-classified with [BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi).

## Quick start (Docker)

```bash
make build
make bash
```

Inside the container, all commands are run from `/usr/src/app/` (the `app/` directory is volume-mounted). The `data/` directory is mounted at `/usr/src/data/`.

## Quick start (local)

Requires [uv](https://docs.astral.sh/uv/), plus the following binaries on `PATH`: `sintax_rs`, `fastq_rs`, `usearch`, `blastn`, `makeblastdb`.

```bash
uv sync
cd app
```

All commands below should be run from the `app/` directory using `uv run`.

## Building the database

Download and format the [EMU](https://github.com/treangenlab/emu) reference database:

```bash
# Docker:
python database.py -o /usr/src/data

# Local:
uv run python database.py -o ../data
```

## Running the pipeline

```bash
# Docker:
python main.py -f /usr/src/data/reads.fastq.gz -d /usr/src/data/database.fasta -o /usr/src/data/output

# Local:
uv run python main.py -f ../data/reads.fastq.gz -d ../data/database.fasta -o ../data/output
```

| Flag | Description | Default |
|------|-------------|---------|
| `-f/--fastq` | Input FASTQ file(s) | required |
| `-d/--database` | Reference database FASTA | required |
| `-o/--outdir` | Output directory | required |
| `--min-len` | Min FASTQ read length  | default: `1200` |
| `--max-len` | Max FASTQ read length  | default: `1700` |
| `--max-error` | Max allowed read error | default: `0.05` (5%) |
| `--cluster-pident` | Percent identity to use during usearch read cluster | default: `0.80` (80%) |
| `-s/--sintax_threshold` | Confidence threshold for taxonomic assignment | default: `0.80` (80%) |
| `-k/--kmer-size` | sintax_rs kmer size | default: `15` |
| `-w/--window-size` | sintax_rs window size | default: `7` |
| `--blast` | Re-classify completely unclassified ASVs with BLAST | default: `false` |
| `--blast-pident` | BLAST percent identity threshold (0.0-1.0) | default: `0.70` (70%) |
| `--blast-qcov` | BLAST query coverage threshold (0.0-1.0) | default: `0.70` (70%) |
| `--threads` | Num threads to use | default: `8` |

## Output files

| File | Description |
|------|-------------|
| `asv.fasta` | Amplicon sequence variants |
| `otutab.tsv` | Read counts per ASV |
| `sintax.tsv` | Raw SINTAX results |
| `parsed.tsv` | Parsed classification results |
| `blast.tsv` | Raw BLAST results (when `--blast` is used) |
| `blast_hits.tsv` | Parsed BLAST hits (when `--blast` is used) |
| `sankey.html` | Taxonomy Sankey diagram |
| `abundance_perc.html` | Relative abundance bar chart |

## Interactive Sankey diagram

```bash
# Docker:
python dash_sankey.py --parsed_tsv /usr/src/data/output/<sample>/parsed.tsv

# Local:
uv run python dash_sankey.py --parsed_tsv ../data/output/<sample>/parsed.tsv
```

Opens a Dash app on port `8000` with a threshold slider for exploring classifications.
