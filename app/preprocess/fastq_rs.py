from pathlib import Path

from common.decorator import with_yaspin
from sh import fastq_rs


@with_yaspin("Running preprocessing...")
def preprocess(
    fastq: Path,
    outdir: Path,
    min_len: int = 1200,
    max_len: int = 1700,
    max_error: float = 0.05,
    threads: int = 2,
) -> Path:
    fasta_out = outdir / "preprocess.fasta"

    filter = fastq_rs.bake(
        "filter",
        "--min-len", min_len,
        "--max-len", max_len,
        "--max-error", max_error,
        "-t", threads,
        "-f",
        _piped=True,
    )
    sort = fastq_rs.bake(
        "sort", "--by", "minimizer", "--reverse", "-t", threads, _piped=True
    )
    to_fa = fastq_rs.bake("fq2-fa", "-t", threads)

    # Filter, sort and convert to fasta.
    to_fa("-o", fasta_out, _in=sort(_in=filter(fastq)))

    assert fasta_out.is_file()
    return fasta_out
