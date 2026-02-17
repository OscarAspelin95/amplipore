import argparse
import logging
import os
import shutil
import tempfile
from dataclasses import dataclass
from pathlib import Path

import pandas as pd

from blast.main import run_blast
from classification.main import classify
from classification.results import get_unclassified_asvs, parse_sintax_tsv
from classification.sintax import run_sintax
from cluster.main import cluster
from common.file import _file, get_file_base, write_subset_fasta
from preprocess.fastq_rs import preprocess

log = logging.getLogger(__name__)

ALLOWED_FASTQ_ENDINGS = (".fastq.gz",)
ALLOWED_FASTA_ENDINGS = (".fasta",)


@dataclass(frozen=True)
class PipelineConfig:
    sintax_threshold: float = 0.80
    kmer_size: int = 15
    window_size: int = 7
    blast: bool = False
    blast_pident: float = 0.70
    blast_qcov: float = 0.70
    cluster_pident: float = 0.80
    threads: int = 8


def run_sample(
    fastq: Path,
    database: Path,
    outdir: Path,
    cfg: PipelineConfig,
) -> tuple[Path, Path, Path, Path | None, list[Path]]:
    # Preprocess and convert fastq to fasta.
    fasta = preprocess(
        fastq,
        outdir,
        threads=cfg.threads,
    )

    # Cluster reads into ASVs.
    asv_fasta, otutab_tsv = cluster(
        fasta,
        outdir,
        pident=cfg.cluster_pident,
        threads=cfg.threads,
    )

    # Clean up intermediary files.
    fasta.unlink(missing_ok=True)
    (outdir / "centroids.fasta").unlink(missing_ok=True)

    # Run SINTAX classification.
    sintax_tsv = run_sintax(asv_fasta, database, outdir, cfg.kmer_size, cfg.window_size)

    # Optionally run BLAST on ASVs that SINTAX couldn't classify.
    blast_df = None
    blast_tsv = None

    if cfg.blast:
        parsed_df = parse_sintax_tsv(sintax_tsv, cfg.sintax_threshold)
        unclassified_asvs = get_unclassified_asvs(parsed_df)

        if unclassified_asvs:
            log.info(f"Running BLAST on {len(unclassified_asvs)} unclassified ASVs.")
            unclassified_fasta = outdir / "unclassified_asvs.fasta"
            write_subset_fasta(asv_fasta, unclassified_asvs, unclassified_fasta)

            blast_df = run_blast(
                unclassified_fasta,
                database,
                outdir,
                pident=cfg.blast_pident,
                qcov=cfg.blast_qcov,
                threads=cfg.threads,
            )
            blast_tsv = outdir / "blast_hits.tsv"
            blast_df.to_csv(blast_tsv, sep="\t", index=False)

            unclassified_fasta.unlink(missing_ok=True)
        else:
            log.info("No unclassified ASVs — skipping BLAST.")

    # Generate results (merges BLAST if available) and visualizations.
    agg_df, plots = classify(
        sintax_tsv, otutab_tsv, cfg.sintax_threshold, outdir, blast_df=blast_df
    )

    agg_tsv = outdir / "agg.tsv"
    agg_df.to_csv(agg_tsv, sep="\t", index=False)

    return asv_fasta, otutab_tsv, agg_tsv, blast_tsv, plots


def copy(file_paths: list[Path | None], target_dir: Path) -> None:
    for f in file_paths:
        if f is not None:
            shutil.copy(f, target_dir)


def main(
    fastq: Path,
    database: Path,
    outdir: Path,
    cfg: PipelineConfig,
) -> None:
    with tempfile.TemporaryDirectory() as temp_dir:
        asv_fasta, otutab_tsv, agg_tsv, blast_tsv, plots = run_sample(
            fastq, database, Path(temp_dir), cfg
        )

        copy([asv_fasta, otutab_tsv, agg_tsv, blast_tsv, *plots], outdir)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="amplipore: 16S amplicon classification pipeline"
    )

    # Required arguments.
    parser.add_argument("-f", "--fastq", help="Path to fastq file", required=True)
    parser.add_argument(
        "-d",
        "--database",
        help="Path to database fasta (default: $DB_PATH)",
        default=None,
    )
    parser.add_argument("-o", "--outdir", required=True)

    # Clustering (usearch).
    parser.add_argument(
        "--cluster-pident",
        type=float,
        default=0.80,
        help="Clustering identity threshold (default: 0.80)",
    )

    # Classification (sintax_rs).
    parser.add_argument(
        "-s",
        "--sintax-threshold",
        type=float,
        default=0.80,
        help="SINTAX confidence threshold (default: 0.80)",
    )
    parser.add_argument(
        "-k",
        "--kmer-size",
        type=int,
        default=15,
        help="SINTAX kmer size (default: 15)",
    )
    parser.add_argument(
        "-w",
        "--window-size",
        type=int,
        default=7,
        help="SINTAX window size (default: 7)",
    )

    # BLAST fallback.
    parser.add_argument(
        "--blast", action="store_true", help="Run BLAST on unclassified ASVs"
    )
    parser.add_argument(
        "--blast-pident",
        type=float,
        default=0.70,
        help="BLAST percent identity threshold (default: 0.70)",
    )
    parser.add_argument(
        "--blast-qcov",
        type=float,
        default=0.70,
        help="BLAST query coverage threshold (default: 0.70)",
    )

    # Shared.
    parser.add_argument(
        "-t", "--threads", type=int, default=8, help="Number of threads (default: 8)"
    )

    args = parser.parse_args()

    fastq = _file(args.fastq, ALLOWED_FASTQ_ENDINGS)
    db_arg = args.database or os.environ.get("DB_PATH")
    if db_arg is None:
        parser.error(
            "--database is required when DB_PATH environment variable is not set"
        )
    database = _file(db_arg, ALLOWED_FASTA_ENDINGS)

    outdir = Path(args.outdir)
    outdir.mkdir(exist_ok=True)

    cfg = PipelineConfig(
        sintax_threshold=args.sintax_threshold,
        kmer_size=args.kmer_size,
        window_size=args.window_size,
        blast=args.blast,
        blast_pident=args.blast_pident,
        blast_qcov=args.blast_qcov,
        cluster_pident=args.cluster_pident,
        threads=args.threads,
    )

    main(fastq, database, outdir, cfg)
