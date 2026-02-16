from pathlib import Path
from sh import usearch


def usearch_cluster(fasta: Path, outdir: Path, pident: float = 0.80, threads: int = 8) -> Path:
    centroid_fasta = outdir / "centroids.fasta"
    usearch(
        "-cluster_fast",
        fasta,
        "-id",
        pident,
        "-centroids",
        centroid_fasta,
        "-strand",
        "both",
        "-sizeout",
        "-threads",
        threads,
    )

    assert centroid_fasta.is_file()
    return centroid_fasta
