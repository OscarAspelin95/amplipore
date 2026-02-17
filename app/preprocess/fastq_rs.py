from pathlib import Path

from common.decorator import with_yaspin
from sh import fastq_rs


@with_yaspin("Running preprocessing...")
def preprocess(
    fastq: Path,
    outdir: Path,
    threads: int = 2,
) -> Path:
    fasta_out = outdir / "preprocess.fasta"

    sort = fastq_rs.bake(
        "sort", "--by", "minimizer", "--reverse", "-t", threads, _piped=True
    )
    to_fa = fastq_rs.bake("fq2-fa", "-t", threads)

    # Sort and convert to fasta.
    to_fa("-o", fasta_out, _in=sort("-f", fastq))

    assert fasta_out.is_file()
    return fasta_out
