from pathlib import Path

import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from plotly.graph_objects import Figure
from common.taxonomy import Levels


def parse_sintax_tsv(tsv: Path, threshold: float) -> pd.DataFrame:
    # ---
    df = pd.read_csv(tsv, sep="\t", names=["asv", "ref", "num_hits", "iteration"])

    # ---
    acc_id_tax = (
        df["ref"]
        .str.split(";", expand=True)
        .rename(columns={0: "accession", 1: "tax_id", 2: "taxonomy"})
        .replace(regex=r"accession=|tax_id=|taxonomy=", value="")
    )
    df = pd.concat([df, acc_id_tax], axis=1).drop(
        columns=["ref", "num_hits", "iteration"]
    )

    # ---
    full_tax_df = (
        df["taxonomy"]
        .str.split("|", expand=True)
        .replace(value="", regex=r"^d:|^p:|^c:|^o:|^f:|^g:|^s:")
        .rename(columns=Levels.pandas_column_rename())
    )
    df = pd.concat([df, full_tax_df], axis=1).drop(columns="taxonomy")

    # ---
    lst = []
    for asv, asv_subset in df.groupby(by="asv"):
        num_iterations = len(asv_subset)

        for level in Levels.as_list():
            (best_hit_at_level, best_absolute_score_at_level) = (
                asv_subset[level]
                .value_counts()
                .reset_index()
                .sort_values(by="count", ascending=False)
                .iloc[0]
            )
            best_relative_score = best_absolute_score_at_level / num_iterations

            lst.append(
                [
                    asv,
                    level,
                    best_hit_at_level
                    if best_relative_score >= threshold
                    else "unclassified",
                    best_relative_score,
                ]
            )

    result_df = pd.DataFrame(lst, columns=["asv", "level", "hit", "score"])

    return result_df


def get_sankey_fig(df: pd.DataFrame) -> Figure:
    levels = Levels.as_list()[::-1]
    edges = []

    for i in range(len(levels) - 1):
        src_col, tgt_col = levels[i], levels[i + 1]

        edges += df[[src_col, tgt_col]].values.tolist()

    # Convert to DataFrame for counting connections
    edges_df = pd.DataFrame(edges, columns=["source", "target"])
    edges_count = edges_df.value_counts().reset_index(name="value")

    # Build unique nodes
    all_nodes = pd.Index(
        pd.concat([edges_count["source"], edges_count["target"]]).unique()
    )

    # Map labels to indices
    node_map = {label: i for i, label in enumerate(all_nodes)}

    # Build Sankey link data
    sources = edges_count["source"].map(node_map)
    targets = edges_count["target"].map(node_map)
    values = edges_count["value"]

    # Create Sankey diagram
    fig = go.Figure(
        go.Sankey(
            node=dict(
                pad=15,
                thickness=20,
                label=all_nodes.tolist(),
            ),
            link=dict(source=sources, target=targets, value=values),
        )
    )

    fig.update_layout(title_text="Taxonomy Sankey Diagram", font_size=10)
    return fig


def pivot_df(df: pd.DataFrame) -> pd.DataFrame:
    pivoted_df = (
        df.pivot(columns="level", index="asv", values="hit").reset_index().dropna()
    )

    for level in Levels.as_list():
        pivoted_df[level] = pivoted_df[level].apply(lambda x: f"{level[0]}:{x}")

    return pivoted_df


def sankey_diagram(parsed_df: pd.DataFrame, outdir: Path):
    pivoted_df = pivot_df(parsed_df)
    fig = get_sankey_fig(pivoted_df)
    fig.write_html(outdir / "sankey.html")


def get_unclassified_asvs(parsed_df: pd.DataFrame) -> set[str]:
    """Return ASV IDs where all 7 taxonomy levels are 'unclassified'."""
    levels = Levels.as_list()
    unclassified = set()

    for asv, group in parsed_df.groupby("asv"):
        hits = group.set_index("level")["hit"]
        if all(hits.get(level, "unclassified") == "unclassified" for level in levels):
            unclassified.add(asv)

    return unclassified


def merge_blast_results(parsed_df: pd.DataFrame, blast_df: pd.DataFrame) -> pd.DataFrame:
    """Replace 'unclassified' entries in SINTAX results with BLAST taxonomy for matched ASVs."""
    if blast_df.empty:
        return parsed_df

    levels = Levels.as_list()
    new_rows = []

    for _, row in blast_df.iterrows():
        asv = row["asv"]
        for level in levels:
            value = row.get(level, "unclassified")
            new_rows.append({"asv": asv, "level": level, "hit": value, "score": None})

    blast_parsed = pd.DataFrame(new_rows)

    # Remove the unclassified SINTAX rows for ASVs that BLAST classified.
    blast_asvs = set(blast_df["asv"])
    kept = parsed_df[~parsed_df["asv"].isin(blast_asvs)]

    return pd.concat([kept, blast_parsed], ignore_index=True)


def get_results(
    parsed_df: pd.DataFrame,
    otutab_tsv: Path,
    outdir: Path,
) -> pd.DataFrame:
    parsed_df.to_csv(outdir / "parsed.tsv", sep="\t", index=False)

    sankey_diagram(parsed_df, outdir)

    otutab_df = pd.read_csv(otutab_tsv, sep="\t")
    total_reads = otutab_df["reads"].sum()

    agg_df = (
        parsed_df.merge(otutab_df, on="asv", how="left")
        .groupby(by=["level", "hit"])["reads"]
        .sum()
        .reset_index()
    )

    agg_df["abundance"] = agg_df["reads"].apply(lambda x: x / total_reads)

    # Sort to ensure proper abundance positioning.
    agg_df = agg_df.sort_values(by=["level", "abundance"], ascending=[False, False])

    agg_df["level"] = pd.Categorical(
        agg_df["level"],
        categories=Levels.as_list(),
        ordered=True,
    )

    fig = px.bar(
        agg_df,
        y="abundance",
        x="level",
        color="hit",
        barmode="relative",
        color_discrete_sequence=px.colors.qualitative.Antique_r,
    )
    fig.write_html(outdir / "abundance_perc.html")

    return agg_df
