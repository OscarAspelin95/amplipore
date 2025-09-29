from common.decorator import with_yaspin
from pathlib import Path
from Bio import SeqIO
import re
import pandas as pd


@with_yaspin("Getting asvs...")
def get_otutab(centroids: Path, outdir: Path) -> tuple[Path, Path]:
    size_pat = re.compile(r"size=(?P<size>\d+);$")

    asv_fasta = outdir / "asv.fasta"
    otutab_tsv = outdir / "otutab.tsv"

    recs, otutab_lst = [], []

    for i, rec in enumerate(SeqIO.parse(centroids, "fasta")):
        if (match := size_pat.search(rec.description)) is None:
            raise ValueError(f"Size missing from {rec.description}")

        size = int(match.groupdict()["size"])
        asv = f"asv_{i}"

        rec.id = asv
        rec.description = ""

        if size <= 1:
            break

        #
        otutab_lst.append([asv, size])
        recs.append(rec)

    with asv_fasta.open("w") as f:
        SeqIO.write(recs, f, "fasta")
    assert asv_fasta.is_file()

    otutab_df = pd.DataFrame(otutab_lst, columns=["asv", "reads"])
    otutab_df.to_csv(otutab_tsv, sep="\t", index=False)
    assert otutab_tsv.is_file()

    return (asv_fasta, otutab_tsv)
