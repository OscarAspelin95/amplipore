from pathlib import Path
from .results import get_results, parse_sintax_tsv, merge_blast_results
import pandas as pd
from common.decorator import with_yaspin


@with_yaspin("Generating classification results...")
def classify(
    sintax_tsv: Path,
    otutab_tsv: Path,
    sintax_threshold: float,
    outdir: Path,
    blast_df: pd.DataFrame | None = None,
) -> pd.DataFrame:
    parsed_df = parse_sintax_tsv(sintax_tsv, sintax_threshold)

    if blast_df is not None:
        parsed_df = merge_blast_results(parsed_df, blast_df)

    return get_results(parsed_df, otutab_tsv, outdir)
