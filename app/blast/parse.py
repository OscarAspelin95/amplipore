import pandas as pd
from pathlib import Path
import logging
from .config import ParseConfig

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


def get_best_hit_per_asv(s: pd.DataFrame) -> pd.DataFrame:
    return s.sort_values(by=["pident", "perc_aln"], ascending=[False, False]).iloc[0]


def parse_blast_tsv(blast_tsv: Path) -> pd.DataFrame:
    cfg = ParseConfig()
    df = pd.read_csv(blast_tsv, sep="\t", names=COLS)

    df = df.assign(
        perc_aln=df["alignment_len"] / df[["subject_len", "query_len"]].min(axis=1),
        # Note, it is rather inefficient to convert to fraction, so for performance reasons
        # we should actually stick with percent.
        pident=df["pident"] / 100,
    )

    # Remove low quality hits.
    df = df.query(f"pident >= {cfg.pident} and perc_aln >= {cfg.perc_aln}")

    # Get the best database hit per asv.
    df = df.groupby(by="subject_id").apply(get_best_hit_per_asv).reset_index(drop=True)
    assert df["subject_id"].is_unique

    # Extract species, remove all other taxonomy.
    df["species"] = df["query_id"].str.extract(r"s:(?P<species>.*)$")
    assert df["species"].isna().sum() == 0

    # NOTE - this log warning probably won't show due to yaspin.
    if df.empty:
        log.warning("No valid BLAST hits found.")

    return df[["subject_id", "species", "pident", "perc_aln"]]
