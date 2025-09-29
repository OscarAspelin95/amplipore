from pydantic import BaseModel
from pathlib import Path
from sh import usearch
from common.decorator import with_yaspin


class UsearchConfig(BaseModel):
    threads: int = 8
    pident: float = 0.80


@with_yaspin("Running read clustering...")
def usearch_cluster(fasta: Path, outdir: Path) -> Path:
    cfg = UsearchConfig()

    centroid_fasta = outdir / "centroids.fasta"
    usearch(
        "-cluster_fast",
        fasta,
        "-id",
        cfg.pident,
        "-centroids",
        centroid_fasta,
        "-strand",
        "both",
        "-sizeout",
        "-threads",
        cfg.threads,
    )

    assert centroid_fasta.is_file()
    return centroid_fasta
