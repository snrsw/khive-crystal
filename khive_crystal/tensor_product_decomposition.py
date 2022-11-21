from copy import deepcopy
from typing import Callable, List, Optional

from khive_crystal.khive import KHive
from khive_crystal.utils import get_length


def decompose(H: List[KHive]) -> KHive:
    pass


def theta(H: List[KHive]) -> List[KHive]:
    pass


def iota(H: KHive) -> Optional[KHive]:
    """iota(H)

    Args:
        H (KHive): A K-hive.

    Returns:
        Optional[KHive]: iota(H)

    Examples:
        >>> H: KHive = KHive(n=3, alpha=[3, 2, 0], beta=[3, 2, 0], gamma=[0, 0, 0], Uij=[[0, 0], [0]])
        >>> iota(H=H)
        KHive(n=3, alpha=[3, 1, 0], beta=[3, 1, 0], gamma=[0, 0, 0], Uij=[[0, 0], [0]])
        >>> H: KHive = KHive(n=3, alpha=[3, 2, 0], beta=[2, 2, 1], gamma=[0, 0, 0], Uij=[[1, 0], [1]])
        >>> iota(H=H)
        KHive(n=3, alpha=[3, 1, 0], beta=[2, 1, 1], gamma=[0, 0, 0], Uij=[[1, 0], [1]])
        >>> H: KHive = KHive(n=3, alpha=[3, 2, 0], beta=[3, 0, 2], gamma=[0, 0, 0], Uij=[[0, 0], [2]])
        >>> iota(H=H)
        KHive(n=3, alpha=[3, 1, 0], beta=[3, 0, 1], gamma=[0, 0, 0], Uij=[[0, 0], [1]])
        >>> H: KHive = KHive(n=3, alpha=[1, 0, 0], beta=[1, 0, 0], gamma=[0, 0, 0], Uij=[[0, 0], [0]])
        >>> iota(H=H)
    """
    return Iota(H=H).run()


class Iota:
    """This class has methods for iota.
    The entry point of this methods is "run". The above function iota is a wrapper of run,
    then use iota instead of using this class directly.
    """
    def __init__(self, H: KHive) -> None:
        self.H: KHive = H
        self.l_alpha: int = get_length(alpha=H.alpha)
        self.l_alpha_as_index: int = get_length(alpha=H.alpha) - 1
        self.j_iota: int = self.compute_j_iota()
        self.j_iota_as_index: int = self.compute_j_iota() - 1

    def compute_j_iota(self) -> int:
        """Compute j_{iota(H)} = min{j in [n] | U_{l(alpha), j} != 0}.

        Returns:
            int: j_{iota(H)}

        Examples:
        >>> H: KHive = KHive(n=3, alpha=[3, 2, 0], beta=[3, 2, 0], gamma=[0, 0, 0], Uij=[[0, 0], [0]])
        >>> Iota(H=H).compute_j_iota()
        2
        >>> H: KHive = KHive(n=3, alpha=[3, 2, 0], beta=[2, 2, 1], gamma=[0, 0, 0], Uij=[[1, 0], [1]])
        >>> Iota(H=H).compute_j_iota()
        2
        >>> H: KHive = KHive(n=3, alpha=[3, 2, 0], beta=[3, 0, 2], gamma=[0, 0, 0], Uij=[[0, 0], [2]])
        >>> Iota(H=H).compute_j_iota()
        3
        >>> H: KHive = KHive(n=3, alpha=[1, 0, 0], beta=[1, 0, 0], gamma=[0, 0, 0], Uij=[[0, 0], [0]])
        >>> Iota(H=H).compute_j_iota()
        1
        """
        j_iota_as_index_in_full_Uij: int = [
            uij == 0
            for uij in self.H.full_Uij[self.l_alpha_as_index]
        ].index(False)
        return j_iota_as_index_in_full_Uij + self.l_alpha

    def compute_alpha(self) -> List[int]:
        """Compute the alpha of iota(H).

        Returns:
            List[int]: the alpha of iota(H).

        Examples:
        >>> H: KHive = KHive(n=3, alpha=[3, 2, 0], beta=[3, 2, 0], gamma=[0, 0, 0], Uij=[[0, 0], [0]])
        >>> Iota(H=H).compute_alpha()
        [3, 1, 0]
        >>> H: KHive = KHive(n=3, alpha=[3, 0, 0], beta=[3, 0, 0], gamma=[0, 0, 0], Uij=[[0, 0], [0]])
        >>> Iota(H=H).compute_alpha()
        [2, 0, 0]
        >>> H: KHive = KHive(n=3, alpha=[1, 0, 0], beta=[1, 0, 0], gamma=[0, 0, 0], Uij=[[0, 0], [0]])
        >>> Iota(H=H).compute_alpha()
        [0, 0, 0]
        """
        return [alpha_i if i != self.l_alpha_as_index else alpha_i - 1 for i, alpha_i in enumerate(self.H.alpha)]

    def compute_beta(self) -> List[int]:
        """Compute the beta of iota(H).

        Returns:
            List[int]: the beta of iota(H).

        Examples:
        >>> H: KHive = KHive(n=3, alpha=[3, 2, 0], beta=[3, 2, 0], gamma=[0, 0, 0], Uij=[[0, 0], [0]])
        >>> Iota(H=H).compute_beta()
        [3, 1, 0]
        >>> H: KHive = KHive(n=3, alpha=[3, 2, 0], beta=[3, 0, 2], gamma=[0, 0, 0], Uij=[[0, 0], [2]])
        >>> Iota(H=H).compute_beta()
        [3, 0, 1]
        >>> H: KHive = KHive(n=3, alpha=[1, 0, 0], beta=[1, 0, 0], gamma=[0, 0, 0], Uij=[[0, 0], [0]])
        >>> Iota(H=H).compute_beta()
        [0, 0, 0]
        """
        return [
            beta_i if i != self.j_iota_as_index else beta_i - 1 for i, beta_i in enumerate(self.H.beta)
        ]

    def compute_Uij(self) -> List[List[int]]:
        """Compute the Uij of iota(H).

        Returns:
            List[List[int]]: the Uij of iota(H).

        Examples:
        >>> H: KHive = KHive(n=3, alpha=[3, 2, 0], beta=[3, 2, 0], gamma=[0, 0, 0], Uij=[[0, 0], [0]])
        >>> Iota(H=H).compute_Uij()
        [[0, 0], [0]]
        >>> H: KHive = KHive(n=3, alpha=[3, 2, 0], beta=[2, 2, 1], gamma=[0, 0, 0], Uij=[[1, 0], [1]])
        >>> Iota(H=H).compute_Uij()
        [[1, 0], [1]]
        >>> H: KHive = KHive(n=3, alpha=[3, 2, 0], beta=[3, 0, 2], gamma=[0, 0, 0], Uij=[[0, 0], [2]])
        >>> Iota(H=H).compute_Uij()
        [[0, 0], [1]]
        >>> H: KHive = KHive(n=3, alpha=[2, 1, 0], beta=[1, 1, 1], gamma=[0, 0, 0], Uij=[[1, 0], [1]])
        >>> Iota(H=H).compute_Uij()
        [[1, 0], [0]]
        """
        if self.l_alpha == self.j_iota:
            return deepcopy(self.H.Uij)
        else:
            return [
                [
                    uij
                    if (i, j) != (self.l_alpha_as_index, self.j_iota_as_index-self.l_alpha) else uij - 1 for j, uij in enumerate(ui)]
                for i, ui in enumerate(self.H.Uij)
            ]

    def run(self) -> Optional[KHive]:
        """Compute iota(H).

        Returns:
            KHive: _description_

        Examples:
        >>> H: KHive = KHive(n=3, alpha=[3, 2, 0], beta=[3, 2, 0], gamma=[0, 0, 0], Uij=[[0, 0], [0]])
        >>> Iota(H=H).run()
        KHive(n=3, alpha=[3, 1, 0], beta=[3, 1, 0], gamma=[0, 0, 0], Uij=[[0, 0], [0]])
        >>> H: KHive = KHive(n=3, alpha=[3, 2, 0], beta=[2, 2, 1], gamma=[0, 0, 0], Uij=[[1, 0], [1]])
        >>> Iota(H=H).run()
        KHive(n=3, alpha=[3, 1, 0], beta=[2, 1, 1], gamma=[0, 0, 0], Uij=[[1, 0], [1]])
        >>> H: KHive = KHive(n=3, alpha=[3, 2, 0], beta=[3, 0, 2], gamma=[0, 0, 0], Uij=[[0, 0], [2]])
        >>> Iota(H=H).run()
        KHive(n=3, alpha=[3, 1, 0], beta=[3, 0, 1], gamma=[0, 0, 0], Uij=[[0, 0], [1]])
        >>> H: KHive = KHive(n=3, alpha=[1, 0, 0], beta=[1, 0, 0], gamma=[0, 0, 0], Uij=[[0, 0], [0]])
        >>> Iota(H=H).run()
        """
        n: int = self.H.n
        alpha: List[int] = self.compute_alpha()
        beta: List[int] = self.compute_beta()
        gamma: List[int] = [0 for _ in range(self.H.n)]
        Uij: List[List[int]] = self.compute_Uij()
        if sum(alpha) == 0:
            return None
        else:
            return KHive(n=n, alpha=alpha, beta=beta, gamma=gamma, Uij=Uij)

def rho(a: int) -> Callable[[KHive], KHive]:
    def _rho(H: KHive) -> KHive:
        pass
    pass
