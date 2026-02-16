import argparse
import logging
from pathlib import Path
from common.file import _file, get_file_base, write_subset_fasta
from preprocess.fastq_rs import preprocess
from cluster.main import cluster
from classification.main import classify
from classification.results import parse_sintax_tsv, get_unclassified_asvs
from classification.sintax import run_sintax
import pandas as pd
from blast.main import run_blast

log = logging.getLogger(__name__)

ALLOWED_FASTQ_ENDINGS = (".fastq.gz",)
ALLOWED_FASTA_ENDINGS = (".fasta",)


def run_sample(
    fastq: Path,
    database: Path,
    sintax_threshold: float,
    blast: bool,
    blast_pident: float,
    blast_qcov: float,
    outdir: Path,
) -> pd.DataFrame:
    sample_name = get_file_base(fastq, ALLOWED_FASTQ_ENDINGS)
    sample_dir = outdir / sample_name
    sample_dir.mkdir(exist_ok=True)

    # Preprocess and convert fastq to fasta.
    fasta = preprocess(fastq, sample_dir)

    # Cluster reads into asvs.
    asv_fasta, otutab_tsv = cluster(fasta, sample_dir)

    # Clean up intermediary files from preprocessing and clustering.
    fasta.unlink(missing_ok=True)
    (sample_dir / "centroids.fasta").unlink(missing_ok=True)

    # Run SINTAX classification.
    sintax_tsv = run_sintax(asv_fasta, database, sample_dir)

    # Optionally run BLAST on ASVs that SINTAX couldn't classify.
    blast_df = None
    if blast:
        parsed_df = parse_sintax_tsv(sintax_tsv, sintax_threshold)
        unclassified_asvs = get_unclassified_asvs(parsed_df)

        if unclassified_asvs:
            log.info(
                f"Running BLAST on {len(unclassified_asvs)} unclassified ASVs."
            )
            unclassified_fasta = sample_dir / "unclassified_asvs.fasta"
            write_subset_fasta(asv_fasta, unclassified_asvs, unclassified_fasta)

            blast_df = run_blast(unclassified_fasta, database, sample_dir, blast_pident, blast_qcov)
            blast_df.to_csv(sample_dir / "blast_hits.tsv", sep="\t", index=False)

            unclassified_fasta.unlink(missing_ok=True)
        else:
            log.info("No unclassified ASVs — skipping BLAST.")

    # Generate results (merges BLAST if available) and visualizations.
    agg_df = classify(sintax_tsv, otutab_tsv, sintax_threshold, sample_dir, blast_df=blast_df)
    agg_df["sample_name"] = sample_name

    return agg_df


def main(
    fastqs: list[Path],
    database: Path,
    sintax_threshold: float,
    blast: bool,
    blast_pident: float,
    blast_qcov: float,
    outdir: Path,
) -> None:
    for fastq in fastqs:
        run_sample(fastq, database, sintax_threshold, blast, blast_pident, blast_qcov, outdir)


if __name__ == "__main__":
    # logging.basicConfig(level=logging.DEBUG, handlers=[logging.StreamHandler()])

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-f", "--fastq", nargs="+", help="Path to fastq file(s).", required=True
    )
    parser.add_argument(
        "-d", "--database", help="Path to database fasta", required=True
    )
    parser.add_argument("-o", "--outdir", required=True)
    parser.add_argument(
        "-s", "--sintax_threshold", type=float, required=False, default=0.80
    )
    parser.add_argument(
        "--blast", action="store_true", help="Run additional classifiction with BLAST"
    )
    parser.add_argument(
        "--blast-pident",
        type=float,
        required=False,
        default=0.70,
        help="BLAST percent identity threshold (0.0-1.0)",
    )
    parser.add_argument(
        "--blast-qcov",
        type=float,
        required=False,
        default=0.70,
        help="BLAST query coverage threshold (0.0-1.0)",
    )
    args = parser.parse_args()

    fastq = [_file(fastq, ALLOWED_FASTQ_ENDINGS) for fastq in args.fastq]
    database = _file(args.database, ALLOWED_FASTA_ENDINGS)
    sintax_threshold = args.sintax_threshold

    # Main output directory.
    outdir = Path(args.outdir)
    outdir.mkdir(exist_ok=True)

    main(
        fastq,
        database,
        sintax_threshold,
        args.blast,
        args.blast_pident,
        args.blast_qcov,
        outdir,
    )
