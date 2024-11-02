import math

import plotly.graph_objects as go

from khive_crystal.khive import KHive


def view(H: KHive | list[KHive]) -> go.Figure:
    """Plot self.H to self.fig. If self.H is a KHive, plot it as origin = (0, 0).
    If self.H is a list of KHive, plot it as tensor products of KHive.
    This is a wrapper of the entry point of the class View.

    Args:
        H (Union[KHive, List[KHive]]):

    Returns:
        go.Figure:
    """
    return View(H=H).run()


class View:
    """This class has methods to view KHive or a list of KHive."""

    def __init__(self, H: KHive | list[KHive]) -> None:
        """init

        Args:
            H (Union[KHive, List[KHive]]): KHive or a list of KHive to plot.
        """
        self.H: KHive | list[KHive] = H
        self.fig: go.Figure = self.setup_figure()
        self.center = self.get_center()

    def get_center(self) -> float:
        """Get center axis of x.

        Returns:
            float: center
        """
        if isinstance(self.H, list):
            return self.H[0].n / 2
        elif isinstance(self.H, KHive):
            return self.H.n / 2

    def setup_figure(self) -> go.Figure:
        """Setup go.Figure

        Returns:
            go.Figure: _description_
        """
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

    def move_origin(self, axis: list[float], origin: int) -> list[float]:
        """Move origin of a graph.

        Args:
            axis (List[float]): Xaxis
            origin (int): Origin

        Returns:
            List[float]: Xaxis moved by origin.
        """
        return [_ + origin for _ in axis]

    def draw_lines_parallel_to_right_edge(self, H: KHive, origin: int) -> None:
        self.fig.add_traces(
            [
                go.Scatter(
                    x=self.move_origin(axis=[_, _ / 2], origin=origin),
                    y=[0, _ / 2],
                    mode="lines",
                    line={"color": "grey"},
                    name=None,
                )
                for _ in range(H.n + 1)
            ]
        )

    def draw_gamma(self, H: KHive, origin: int) -> None:
        self.fig.add_traces(
            [
                go.Scatter(
                    x=self.move_origin(axis=[(_ + 0.2) / 2], origin=origin),
                    y=[(_ + 0.7) / 2],
                    text=[H.gamma[_]],
                    mode="text",
                    textfont={"color": "black"},
                    name=None,
                )
                for _ in range(H.n)
            ]
        )

    def draw_lines_parallel_to_left_edge(self, H: KHive, origin: int) -> None:
        self.fig.add_traces(
            [
                go.Scatter(
                    x=self.move_origin(axis=[_, self.center + _ * 0.5], origin=origin),
                    y=[0, self.center - _ * 0.5],
                    mode="lines",
                    line={"color": "grey"},
                    name=None,
                )
                for _ in range(H.n)
            ]
        )

    def draw_alpha(self, H: KHive, origin: int) -> None:
        self.fig.add_traces(
            [
                go.Scatter(
                    x=self.move_origin(
                        axis=[self.center + (_ + 0.8) * 0.5], origin=origin
                    ),
                    y=[self.center - (_ + 0.3) * 0.5],
                    text=[H.alpha[_]],
                    mode="text",
                    textfont={"color": "black"},
                    name=None,
                )
                for _ in range(H.n)
            ]
        )

    def draw_below_edge(self, H: KHive, origin: int) -> None:
        self.fig.add_traces(
            [
                go.Scatter(
                    x=self.move_origin(axis=[0, H.n], origin=origin),
                    y=[0, 0],
                    mode="lines",
                    line={"color": "grey"},
                    name=None,
                )
            ]
        )

    def draw_beta(self, H: KHive, origin: int) -> None:
        self.fig.add_traces(
            [
                go.Scatter(
                    x=self.move_origin(axis=[0.5 + _], origin=origin),
                    y=[-0.15],
                    text=[H.beta[_]],
                    mode="text",
                    textfont={"color": "black"},
                    name=None,
                )
                for _ in range(H.n)
            ]
        )

    def draw_Uij(self, H: KHive, origin: int) -> None:
        self.fig.add_traces(
            [
                go.Scatter(
                    x=self.move_origin(axis=[i + 1 + 0.5 * j], origin=origin),
                    y=[0.5 + 0.5 * j],
                    text=[Uij],
                    mode="text",
                    textfont={"color": "black"},
                    name=None,
                )
                for i, Ui in enumerate(H.Uij)
                for j, Uij in enumerate(Ui)
            ]
        )

    def plot_otimes(self, origin: int) -> None:
        """Plot otimes at (elf.center * 2 + 0.5 + origin, self.center / 2).
        This function is used when H is a list of Khives.

        Args:
            origin (int): Origin
        """
        self.fig.add_annotation(
            x=self.center * 2 + 0.5 + origin,
            y=self.center / 2,
            text=r"$\otimes$",
            font={"size": self.center * 40},
            showarrow=False,
        )

    def plot(self, H: KHive, origin: int) -> None:
        """Plot a KHive to self.fig.

        Args:
            H (KHive): KHive
            origin (int): origin
        """
        self.draw_below_edge(H=H, origin=origin)
        self.draw_lines_parallel_to_left_edge(H=H, origin=origin)
        self.draw_lines_parallel_to_right_edge(H=H, origin=origin)
        self.draw_alpha(H=H, origin=origin)
        self.draw_beta(H=H, origin=origin)
        self.draw_gamma(H=H, origin=origin)
        self.draw_Uij(H=H, origin=origin)

    def run(self) -> go.Figure:
        """Plot self.H to self.fig. If self.H is a KHive, plot it as origin = (0, 0).
        If self.H is a list of KHive, plot it as tensor products of KHive.

        Returns:
            go.Figure: self.fig
        """
        if isinstance(self.H, list) and len(self.H) == 1:
            self.plot(H=self.H[0], origin=0)
            return self.fig

        if isinstance(self.H, KHive):
            self.plot(H=self.H, origin=0)
            return self.fig

        n: int = self.H[0].n
        for i, Hi in enumerate(self.H):
            origin: int = i + i * n
            self.plot(H=Hi, origin=origin)
            if i != len(self.H) - 1:
                self.plot_otimes(origin=origin)
        return self.fig
