import math

import plotly.graph_objects as go
from khive_crystal.khive import KHive


def view(H: KHive) -> go.Figure:
    return View(H=H).run()


class View:
    def __init__(self, H: KHive) -> None:
        self.H = H
        self.center: float = H.n / 2
        self.fig: go.Figure = self.setup_figure()

    def setup_figure(self) -> go.Figure:
        fig: go.Figure = go.Figure()
        fig.update_layout(
            template="plotly_white",
            xaxis={"visible": False, "showticklabels": False},
            yaxis={
                "visible": False,
                "showticklabels": False,
                "scaleanchor": "x",
                "scaleratio": math.sqrt(3),
            },
            showlegend=False,
            autosize=True,
        )
        return fig

    def draw_lines_parallel_to_right_edge(self) -> None:
        self.fig.add_traces(
            [
                go.Scatter(
                    x=[_, _ / 2],
                    y=[0, _ / 2],
                    mode="lines",
                    line={"color": "grey"},
                    name=None,
                )
                for _ in range(self.H.n + 1)
            ]
        )

    def draw_lines_parallel_to_left_edge(self) -> None:
        self.fig.add_traces(
            [
                go.Scatter(
                    x=[_, self.center + _ * 0.5],
                    y=[0, self.center - _ * 0.5],
                    mode="lines",
                    line={"color": "grey"},
                    name=None,
                )
                for _ in range(self.H.n)
            ]
        )

    def draw_below_edge(self) -> None:
        self.fig.add_traces(
            [
                go.Scatter(
                    x=[0, self.H.n],
                    y=[0, 0],
                    mode="lines",
                    line={"color": "grey"},
                    name=None,
                )
            ]
        )

    def draw_Uij(self) -> None:
        self.fig.add_traces(
            [
                go.Scatter(
                    x=[i + 1 + 0.5 * j],
                    y=[0.5 + 0.5 * j],
                    text=[Uij],
                    mode="text",
                    textfont={"color": "black"},
                    name=None,
                )
                for i, Ui in enumerate(self.H.Uij)
                for j, Uij in enumerate(Ui)
            ]
        )

    def run(self) -> go.Figure:
        self.draw_below_edge()
        self.draw_lines_parallel_to_left_edge()
        self.draw_lines_parallel_to_right_edge()
        self.draw_Uij()
        return self.fig
