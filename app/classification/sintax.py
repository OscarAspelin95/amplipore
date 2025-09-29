from sh import sintax_rs
from pathlib import Path
from pydantic import BaseModel
from common.decorator import with_yaspin


class SintaxConfig(BaseModel):
    num_query_hashes: int = 32
    num_bootstrap: int = 100
    kmer_size: int = 15


@with_yaspin("Running classification...")
def run_sintax(asv_fasta: Path, database: Path, outdir: Path) -> Path:
    cfg = SintaxConfig()

    sintax_tsv = outdir / "sintax.tsv"

    sintax_rs(
        "--query",
        asv_fasta,
        "--database",
        database,
        "--outfile",
        sintax_tsv,
        "--bootstraps",
        cfg.num_bootstrap,
        "--query-hashes",
        cfg.num_query_hashes,
        "--kmer-size",
        cfg.kmer_size,
    )

    assert sintax_tsv.is_file()
    return sintax_tsv
