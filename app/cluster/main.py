from pathlib import Path
from .usearch import usearch_cluster
from .otutab import get_otutab
from common.decorator import with_yaspin


@with_yaspin(progress_text="Generating ASVs...")
def cluster(fasta: Path, outdir: Path) -> tuple[Path, Path]:
    centroids = usearch_cluster(fasta, outdir)

    asv_fasta, otutab_tsv = get_otutab(centroids, outdir)
    return asv_fasta, otutab_tsv
