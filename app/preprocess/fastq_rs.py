from sh import fastq_rs
from common.decorator import with_yaspin
from pydantic import BaseModel
from pathlib import Path


class FastqConfig(BaseModel):
    min_len: int = 1200
    max_len: int = 1700
    max_error: int = 0.05
    threads: int = 2


@with_yaspin("Running preprocessing...")
def preprocess(fastq: Path, outdir: Path):
    cfg = FastqConfig()
    fasta_out = outdir / "preprocess.fasta"

    filter = fastq_rs.bake(
        "filter",
        "--min-len",
        cfg.min_len,
        "--max-len",
        cfg.max_len,
        "--max-error",
        cfg.max_error,
        "-t",
        cfg.threads,
        "-f",
        _piped=True,
    )
    sort = fastq_rs.bake(
        "sort", "--by", "minimizer", "--reverse", "-t", cfg.threads, _piped=True
    )
    to_fa = fastq_rs.bake("fq2-fa", "-t", cfg.threads)

    # Filter, sort and convert to fasta.
    to_fa("-o", fasta_out, _in=sort(_in=filter(fastq)))

    assert fasta_out.is_file()
    return fasta_out
