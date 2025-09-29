import argparse
from pathlib import Path
from Bio import SeqIO
import pandas as pd
import re
from sh import curl
from common.decorator import with_yaspin

# Path to EMU github
EMU_GITHUB = Path(
    "https://raw.githubusercontent.com/treangenlab/emu/master/emu_database/"
)

tax_id_pat = re.compile(r"(?P<tax_id>(\d)+):emu_db:\d+")
accession_pat = re.compile(r"[A-Z]{2}_\d+.\d")


def get_tax_id(rec_id: str) -> int:
    match tax_id_pat.match(rec_id):
        case None:
            raise ValueError(rec_id)

        case _ as tax_id_match:
            tax_id = tax_id_match.groupdict()["tax_id"]
            return tax_id


def get_accession(rec_desc: str) -> str:
    match accession_pat.search(rec_desc):
        case None:
            raise ValueError(rec_desc)
        case _ as accession_match:
            return accession_match.group(0)


def sanitize(s: str) -> str:
    return s.lower().replace(" ", "_")


def get_taxonomy(tax_id: int, df: pd.DataFrame) -> str:
    query = df.query(f"tax_id == {tax_id}")

    if len(query) != 1:
        raise ValueError(f"{tax_id}, {query}")

    query = query.iloc[0]

    _, species, genus, family, order, clas, phylum, _, domain, *_ = query

    taxonomy = (
        f"d:{domain}|p:{phylum}|c:{clas}|o:{order}|f:{family}|g:{genus}|s:{species}"
    )
    sanitized_taxonomy = sanitize(taxonomy)

    assert sanitized_taxonomy.count("|") == 6
    return sanitized_taxonomy


def fetch_fasta(outdir: Path) -> Path:
    emu_fasta = outdir / "emu.fasta"

    curl(EMU_GITHUB / "species_taxid.fasta", "-o", emu_fasta)
    assert emu_fasta.is_file()

    return emu_fasta


def fetch_taxonomy(outdir: Path) -> Path:
    emu_taxonomy = outdir / "emu.tsv"

    curl(EMU_GITHUB / "taxonomy.tsv", "-o", emu_taxonomy)
    assert emu_taxonomy.is_file()

    return emu_taxonomy


@with_yaspin("--- Downloading files...")
def fetch_db(outdir: Path) -> tuple[Path, Path]:
    emu_fasta = fetch_fasta(outdir)
    emu_taxonomy = fetch_taxonomy(outdir)

    return emu_fasta, emu_taxonomy


@with_yaspin("--- Parsing taxonomy...")
def write_db(outdir: Path, emu_fasta: Path, emu_taxonomy: Path) -> Path:
    db_fasta = outdir / "db.fasta"

    df = pd.read_csv(emu_taxonomy, sep="\t")

    lst = []
    with emu_fasta.open("r") as f:
        for rec in SeqIO.parse(f, "fasta"):
            tax_id = get_tax_id(rec.id)
            accession = get_accession(rec.description)
            taxonomy = get_taxonomy(tax_id, df)

            rec.id = f"accession={accession};tax_id={tax_id};taxonomy={taxonomy}"
            rec.description = ""

            lst.append(rec)

    with db_fasta.open("w") as f:
        SeqIO.write(lst, f, "fasta")

    assert db_fasta.is_file()
    return db_fasta


def main(outdir: Path) -> Path:
    emu_fasta, emu_taxonomy = fetch_db(outdir)

    db_fasta = write_db(outdir, emu_fasta, emu_taxonomy)

    emu_fasta.unlink()
    emu_taxonomy.unlink()

    return db_fasta


if __name__ == "__main__":
    # Arguments.
    parser = argparse.ArgumentParser()
    parser.add_argument("-o", "--outdir", help="Output directory", required=True)
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(exist_ok=True, parents=True)

    _ = main(outdir)
