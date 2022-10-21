from dataclasses import dataclass, field

import graphviz

from khive_crystal.crystal_structure import e, epsilon, f, phi
from khive_crystal.khive import KHive
from khive_crystal.khives import KHives


@dataclass(frozen=True)
class CrystalGraph:
    khives: KHives
    G: graphviz.Digraph = field(init=False, default=graphviz.Digraph(format="png"))

    def lower_graph(self, H: KHive) -> None:
        self.G.node(str(H))
        for i in range(1, H.n):
            if phi(i=i)(H) > 0:
                K: KHive = f(i=i)(H)  # type: ignore
                self.G.edge(tail_name=str(H), head_name=str(K), label=i)
                self.lower_graph(H=K)

    def raiging_graph(self, H: KHive) -> None:
        self.G.node(H)
        for i in range(1, H.n):
            if epsilon(i=i)(H) > 0:
                K: KHive = e(i=i)(H)  # type: ignore
                self.G.edge(tail_name=str(K), head_name=str(H), label=i)
                self.raiging_graph(H=K)

    def run(self) -> graphviz.Digraph:
        H: KHive = self.khives.highest_weight_vector()
        self.lower_graph(H=H)
        self.raiging_graph(H=H)
        return self.G


def crystal_graph(khives: KHives) -> graphviz.Digraph:
    crystal_graph = CrystalGraph(khives=khives).run()
    return crystal_graph
