from pathlib import Path
from .usearch import usearch_cluster
from .otutab import get_otutab
from common.decorator import with_yaspin


@with_yaspin(progress_text="Generating ASVs...")
def cluster(fasta: Path, outdir: Path, pident: float = 0.80, threads: int = 8) -> tuple[Path, Path]:
    centroids = usearch_cluster(fasta, outdir, pident=pident, threads=threads)

    asv_fasta, otutab_tsv = get_otutab(centroids, outdir)
    return asv_fasta, otutab_tsv
