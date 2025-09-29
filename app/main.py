import argparse
import logging
from pathlib import Path
from common.file import _file, get_file_base
from preprocess.fastq_rs import preprocess
from cluster.main import cluster
from classification.main import classify
import sys
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parent
sys.path.append(str(PROJECT_ROOT))

log = logging.getLogger(__name__)

ALLOWED_FASTQ_ENDINGS = (".fastq.gz",)
ALLOWED_FASTA_ENDINGS = (".fasta",)


def run_sample(
    fastq: Path,
    database: Path,
    sintax_threshold: float,
    outdir: Path,
) -> pd.DataFrame:
    sample_name = get_file_base(fastq, ALLOWED_FASTQ_ENDINGS)
    sample_dir = outdir / sample_name
    sample_dir.mkdir(exist_ok=True)

    # Preprocess and convert fastq to fasta.
    fasta = preprocess(fastq, sample_dir)

    # Cluster reads into asvs.
    asv_fasta, otutab_tsv = cluster(fasta, sample_dir)

    # Classify asvs.
    agg_df = classify(asv_fasta, otutab_tsv, database, sintax_threshold, sample_dir)
    agg_df["sample_name"] = sample_name

    return agg_df


def main(
    fastqs: list[Path],
    database: Path,
    sintax_threshold: float,
    outdir: Path,
) -> None:
    for fastq in fastqs:
        run_sample(fastq, database, sintax_threshold, outdir)


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
    args = parser.parse_args()

    fastq = [_file(fastq, ALLOWED_FASTQ_ENDINGS) for fastq in args.fastq]
    database = _file(args.database, ALLOWED_FASTA_ENDINGS)
    sintax_threshold = args.sintax_threshold

    # Main output directory.
    outdir = Path(args.outdir)
    outdir.mkdir(exist_ok=True)

    main(fastq, database, sintax_threshold, outdir)
