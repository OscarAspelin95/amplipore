import pandas as pd
from pathlib import Path
import logging
from common.taxonomy import Levels, TAXONOMY_PREFIXES

log = logging.getLogger(__name__)

COLS = [
    "query_id",
    "query_len",
    "query_start",
    "query_end",
    "query_frame",
    "subject_id",
    "subject_len",
    "subject_start",
    "subject_end",
    "subject_frame",
    "alignment_len",
    "pident",
]

TAXONOMY_LEVELS = Levels.as_list()


def get_best_hit_per_asv(s: pd.DataFrame) -> pd.DataFrame:
    return s.sort_values(by=["pident", "perc_aln"], ascending=[False, False]).iloc[0]


def extract_taxonomy(subject_id: str) -> dict[str, str]:
    """Extract taxonomy from reference header format: accession=...;tax_id=...;taxonomy=d:...|p:...|...|s:..."""
    taxonomy = {level: "unclassified" for level in TAXONOMY_LEVELS}

    parts = subject_id.split(";")
    for part in parts:
        if part.startswith("taxonomy="):
            tax_str = part.removeprefix("taxonomy=")
            for field in tax_str.split("|"):
                if ":" in field:
                    prefix, value = field.split(":", 1)
                    if prefix in TAXONOMY_PREFIXES and value:
                        taxonomy[TAXONOMY_PREFIXES[prefix]] = value

    return taxonomy


def parse_blast_tsv(blast_tsv: Path, pident: float = 0.70, perc_aln: float = 0.70) -> pd.DataFrame:
    df = pd.read_csv(blast_tsv, sep="\t", names=COLS)

    if df.empty:
        log.warning("No BLAST hits found.")
        return pd.DataFrame(columns=["asv"] + TAXONOMY_LEVELS + ["pident", "perc_aln"])

    df = df.assign(
        perc_aln=df["alignment_len"] / df[["subject_len", "query_len"]].min(axis=1),
        pident=df["pident"] / 100,
    )

    # Remove low quality hits.
    df = df.query(f"pident >= {pident} and perc_aln >= {perc_aln}")

    if df.empty:
        log.warning("No valid BLAST hits after filtering.")
        return pd.DataFrame(columns=["asv"] + TAXONOMY_LEVELS + ["pident", "perc_aln"])

    # Get the best database hit per ASV (query_id = ASV after fixing direction).
    df = df.groupby(by="query_id").apply(get_best_hit_per_asv).reset_index(drop=True)
    assert df["query_id"].is_unique

    # Extract full taxonomy from subject_id (reference header).
    tax_rows = df["subject_id"].apply(extract_taxonomy).apply(pd.Series)
    df = pd.concat([df, tax_rows], axis=1)

    df = df.rename(columns={"query_id": "asv"})

    return df[["asv"] + TAXONOMY_LEVELS + ["pident", "perc_aln"]]
