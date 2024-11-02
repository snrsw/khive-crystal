import graphviz

from khive_crystal.crystal_structure import e, epsilon, f, phi
from khive_crystal.khive import KHive
from khive_crystal.khives import KHives


class CrystalGraph:
    """This class has methods to generate the crystal graph for a KHives.
    The entry point of this methods is "run". The function crystal_graph is a wrapper of run,
    then use crystal_graph instead of using this class direct.
    """

    def __init__(self, khives: KHives):
        self.khives: KHives = khives
        self.G: graphviz.Digraph = graphviz.Digraph(strict=True)

    def lower_graph(self, H: KHive) -> None:
        self.G.node(str(H))
        for i in range(1, H.n):
            if phi(i=i)(H) > 0:
                K: KHive = f(i=i)(H)  # type: ignore
                self.G.edge(tail_name=str(H), head_name=str(K), label=str(i))
                self.lower_graph(H=K)

    def raiging_graph(self, H: KHive) -> None:
        self.G.node(str(H))
        for i in range(1, H.n):
            if epsilon(i=i)(H) > 0:
                K: KHive = e(i=i)(H)  # type: ignore
                self.G.edge(tail_name=str(K), head_name=str(H), label=str(i))
                self.raiging_graph(H=K)

    def run(self) -> graphviz.Digraph:
        H: KHive = self.khives.highest_weight_vector()
        self.lower_graph(H=H)
        self.raiging_graph(H=H)
        return self.G


def crystal_graph(khives: KHives) -> graphviz.Digraph:
    """Generate the crystal graph for given KHives

    Args:
        khives (KHives): khives

    Returns:
        graphviz.Digraph: the crystal graph for khives
    """
    crystal_graph: graphviz.Digraph = CrystalGraph(khives=khives).run()
    return crystal_graph
