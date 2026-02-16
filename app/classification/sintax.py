from pathlib import Path

from common.decorator import with_yaspin
from sh import sintax_rs


@with_yaspin("Running SINTAX classification...")
def run_sintax(
    asv_fasta: Path,
    database: Path,
    outdir: Path,
    kmer_size: int,
    window_size: int,
    num_bootstrap: int = 100,
    num_query_hashes: int = 32,
) -> Path:
    sintax_tsv = outdir / "sintax.tsv"

    sintax_rs(
        "--query",
        asv_fasta,
        "--database",
        database,
        "--outfile",
        sintax_tsv,
        "--bootstraps",
        num_bootstrap,
        "--query-hashes",
        num_query_hashes,
        "--kmer-size",
        kmer_size,
        "--window-size",
        window_size,
    )

    assert sintax_tsv.is_file()
    return sintax_tsv
