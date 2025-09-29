from pathlib import Path
from .sintax import run_sintax
from .results import get_results
import pandas as pd


def classify(
    asv_fasta: Path,
    otutab_tsv: Path,
    database: Path,
    sintax_threshold: float,
    outdir: Path,
) -> pd.DataFrame:
    sintax_tsv = run_sintax(asv_fasta, database, outdir)

    agg_df = get_results(sintax_tsv, otutab_tsv, sintax_threshold, outdir)

    return agg_df
