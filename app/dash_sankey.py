from pathlib import Path
import pandas as pd
from dash import Dash, dcc, html, Input, Output
import argparse
import os
from classification.results import get_sankey_fig, pivot_df
from common.file import _file


def dash_wrapper(tsv: Path):
    df = pd.read_csv(tsv, sep="\t")

    app = Dash(__name__)

    app.layout = html.Div(
        [
            html.H4("Result Chart"),
            dcc.Graph(id="graph"),
            html.P("Threshold"),
            dcc.Slider(
                id="slider",
                min=round(df["score"].min(), 1),
                max=1.0,
                value=round(df["score"].min(), 1),
                step=0.1,
            ),
        ]
    )

    @app.callback(Output("graph", "figure"), Input("slider", "value"))
    def sankey_diagram(threshold: float):
        local_df = df.copy()

        # Update the dataframe with unclassified when threshold is updated.
        local_df["hit"] = local_df.apply(
            lambda x: x["hit"] if x["score"] >= threshold else "unclassified",
            axis=1,
        )

        local_df = pivot_df(local_df)
        fig = get_sankey_fig(local_df)

        return fig

    app.run(debug=True, port=os.environ.get("DASH_PORT"), host="0.0.0.0")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--parsed_tsv",
        required=True,
        help="Path to the 'parsed.tsv' file",
    )

    args = parser.parse_args()

    tsv = _file(args.parsed_tsv, allowed_file_endings=(".tsv"))

    dash_wrapper(tsv)
