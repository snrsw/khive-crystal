from copy import deepcopy
from typing import Callable, List, Optional, Tuple

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
    """  # noqa: B950
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
        """  # noqa: B950
        j_iota_as_index_in_full_Uij: int = [
            uij == 0 for uij in self.H.full_Uij[self.l_alpha_as_index]
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
        """  # noqa: B950
        return [
            alpha_i if i != self.l_alpha_as_index else alpha_i - 1
            for i, alpha_i in enumerate(self.H.alpha)
        ]

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
        """  # noqa: B950
        return [
            beta_i if i != self.j_iota_as_index else beta_i - 1
            for i, beta_i in enumerate(self.H.beta)
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
        """  # noqa: B950
        if self.l_alpha == self.j_iota:
            return deepcopy(self.H.Uij)
        else:
            return [
                [
                    uij
                    if (i, j)
                    != (self.l_alpha_as_index, self.j_iota_as_index - self.l_alpha)
                    else uij - 1
                    for j, uij in enumerate(ui)
                ]
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
        """  # noqa: B950
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
    """Compute rho_i(H).

    Returns:
        KHive: rho_i(H)

    Examples:
    >>> H: KHive = KHive(n=3, alpha=[3, 2, 0], beta=[3, 2, 0], gamma=[0, 0, 0], Uij=[[0, 0], [0]])
    >>> rho(a=1)(H)
    KHive(n=3, alpha=[4, 2, 0], beta=[4, 2, 0], gamma=[0, 0, 0], Uij=[[0, 0], [0]])
    >>> rho(a=2)(H)
    KHive(n=3, alpha=[4, 2, 0], beta=[3, 3, 0], gamma=[0, 0, 0], Uij=[[1, 0], [0]])
    >>> H: KHive = KHive(n=3, alpha=[3, 2, 0], beta=[2, 2, 1], gamma=[0, 0, 0], Uij=[[0, 1], [0]])
    >>> rho(a=2)(H)
    KHive(n=3, alpha=[3, 3, 0], beta=[2, 3, 1], gamma=[0, 0, 0], Uij=[[1, 0], [1]])
    >>> H: KHive = KHive(n=3, alpha=[3, 2, 0], beta=[1, 2, 2], gamma=[0, 0, 0], Uij=[[2, 0], [2]])
    >>> rho(a=1)(H)
    KHive(n=3, alpha=[2, 1, 0], beta=[1, 1, 1], gamma=[0, 0, 0], Uij=[[1, 0], [1]])
    """  # noqa: B950

    def _rho(H: KHive) -> KHive:
        return Rho(H=H, a=a).run()

    return _rho


class Rho:
    """This class has methods for rho.
    The entry point of this methods is "run". The above function rho is a wrapper of run,
    then use rho instead of using this class directly.
    """

    def __init__(self, H: KHive, a: int) -> None:
        self.H: KHive = H
        self.a: int = a
        self.a_as_index: int = a - 1
        self.path: List[Tuple[int, int]] = self.get_path()

    def get_path(self) -> List[Tuple[int, int]]:
        """Get a path on H, which is used to compute rho(H).

        Returns:
            List[Tuple[int, int]]: Path

        Examples:
            >>> H: KHive = KHive(
            ...     n=4,
            ...     alpha=[6, 4, 1, 0],
            ...     beta=[3, 4, 2, 2],
            ...     gamma=[0, 0, 0, 0],
            ...     Uij=[
            ...         [2, 0, 1],
            ...         [1, 1],
            ...         [0]
            ...     ]
            ... )
            >>> Rho(H=H, a=1).get_path()
            [(0, 1), (1, 1), (1, 2), (2, 2), (2, 3), (3, 3), (3, 5)]
            >>> H: KHive = KHive(
            ...     n=4,
            ...     alpha=[6, 4, 1, 0],
            ...     beta=[3, 4, 2, 2],
            ...     gamma=[0, 0, 0, 0],
            ...     Uij=[
            ...         [2, 0, 1],
            ...         [1, 1],
            ...         [0]
            ...     ]
            ... )
            >>> Rho(H=H, a=3).get_path()
            [(0, 3), (1, 3), (1, 4), (2, 4), (2, 5)]
        """
        path: List[Tuple[int, int]] = []
        path.append((0, self.a))

        while True:
            path.append((path[-1][0] + 1, path[-1][1]))  # even case
            next_j_candidate: List[int] = [
                j
                for j in range(path[-1][1] + 1, self.H.n + 1)
                if self.H.Uij[path[-1][0] - 1][j - path[-1][0] - 1] > 0
            ]

            if next_j_candidate == []:
                path.append((path[-1][0], self.H.n + 1))  # odd case
                break

            path.append((path[-1][0], min(next_j_candidate)))  # odd case

        return path

    def compute_alpha(self) -> List[int]:
        """Compute the alpha of rho_i(H).

        Returns:
            List[int]: the alpha of rho_i(H).

        Examples:
        >>> H: KHive = KHive(n=3, alpha=[3, 2, 0], beta=[3, 2, 0], gamma=[0, 0, 0], Uij=[[0, 0], [0]])
        >>> Rho(H=H, a=1).compute_alpha()
        [4, 2, 0]
        >>> Rho(H=H, a=2).compute_alpha()
        [4, 2, 0]
        >>> H: KHive = KHive(n=3, alpha=[3, 2, 0], beta=[2, 2, 1], gamma=[0, 0, 0], Uij=[[0, 1], [0]])
        >>> Rho(H=H, a=2).compute_alpha()
        [3, 3, 0]
        >>> H: KHive = KHive(n=3, alpha=[3, 2, 0], beta=[1, 2, 2], gamma=[0, 0, 0], Uij=[[2, 0], [2]])
        >>> Rho(H=H, a=1).compute_alpha()
        [2, 1, 0]
        """  # noqa: B950
        i_N: int = self.path[-1][0]
        i_N_as_index: int = i_N - 1
        n_as_index: int = self.H.n - 1

        if i_N < self.H.n:
            return [
                alpha_i + 1 if i == i_N_as_index else alpha_i
                for i, alpha_i in enumerate(self.H.alpha)
            ]
        else:
            return [
                alpha_i - 1 if i != n_as_index else alpha_i
                for i, alpha_i in enumerate(self.H.alpha)
            ]

    def compute_beta(self) -> List[int]:
        """Compute the beta of rho_i(H).

        Returns:
            List[int]: the beta of rho_i(H).

        Examples:
        >>> H: KHive = KHive(n=3, alpha=[3, 2, 0], beta=[3, 2, 0], gamma=[0, 0, 0], Uij=[[0, 0], [0]])
        >>> Rho(H=H, a=1).compute_beta()
        [4, 2, 0]
        >>> Rho(H=H, a=2).compute_beta()
        [3, 3, 0]
        >>> Rho(H=H, a=3).compute_beta()
        [3, 2, 1]
        """  # noqa: B950
        j_0_as_index: int = self.path[0][1] - 1
        i_N: int = self.path[-1][0]

        if i_N < self.H.n:
            return [
                beta_i if i != j_0_as_index else beta_i + 1
                for i, beta_i in enumerate(self.H.beta)
            ]
        else:
            return [
                beta_i - 1 if i != j_0_as_index else beta_i
                for i, beta_i in enumerate(self.H.beta)
            ]

    def compute_Uij(self) -> List[List[int]]:
        """Compute the Uij of rho_i(H).

        Returns:
            List[List[int]]: the Uij of rho_i(H).

        Examples:
        >>> H: KHive = KHive(n=3, alpha=[3, 2, 0], beta=[3, 2, 0], gamma=[0, 0, 0], Uij=[[0, 0], [0]])
        >>> Rho(H=H, a=1).compute_Uij()
        [[0, 0], [0]]
        >>> Rho(H=H, a=2).compute_Uij()
        [[1, 0], [0]]
        >>> H: KHive = KHive(n=3, alpha=[3, 2, 0], beta=[2, 2, 1], gamma=[0, 0, 0], Uij=[[0, 1], [0]])
        >>> Rho(H=H, a=2).compute_Uij()
        [[1, 0], [1]]
        >>> H: KHive = KHive(n=3, alpha=[3, 2, 0], beta=[1, 2, 2], gamma=[0, 0, 0], Uij=[[2, 0], [2]])
        >>> Rho(H=H, a=1).compute_Uij()
        [[1, 0], [1]]
        """  # noqa: B950
        Uij: List[List[int]] = deepcopy(self.H.Uij)

        is_point_on_edge: Callable[[int], bool] = lambda m: (m == 0) or (
            m == len(self.path) - 1
        )
        is_point_on_uii: Callable[[int, int], bool] = lambda i, j: i == j

        for m, p_m in enumerate(self.path):
            i_m_as_index: int = p_m[0] - 1
            j_m_as_index: int = p_m[1] - p_m[0] - 1

            if is_point_on_edge(m) or is_point_on_uii(p_m[0], p_m[1]):
                continue
            elif m % 2 == 0:
                Uij[i_m_as_index][j_m_as_index] -= 1
            else:
                Uij[i_m_as_index][j_m_as_index] += 1

        return Uij

    def run(self) -> KHive:
        """Compute rho_i(H).

        Returns:
            KHive: rho_i(H)

        Examples:
        >>> H: KHive = KHive(n=3, alpha=[3, 2, 0], beta=[3, 2, 0], gamma=[0, 0, 0], Uij=[[0, 0], [0]])
        >>> Rho(H=H, a=1).run()
        KHive(n=3, alpha=[4, 2, 0], beta=[4, 2, 0], gamma=[0, 0, 0], Uij=[[0, 0], [0]])
        >>> Rho(H=H, a=2).run()
        KHive(n=3, alpha=[4, 2, 0], beta=[3, 3, 0], gamma=[0, 0, 0], Uij=[[1, 0], [0]])
        >>> H: KHive = KHive(n=3, alpha=[3, 2, 0], beta=[2, 2, 1], gamma=[0, 0, 0], Uij=[[0, 1], [0]])
        >>> Rho(H=H, a=2).run()
        KHive(n=3, alpha=[3, 3, 0], beta=[2, 3, 1], gamma=[0, 0, 0], Uij=[[1, 0], [1]])
        >>> H: KHive = KHive(n=3, alpha=[3, 2, 0], beta=[1, 2, 2], gamma=[0, 0, 0], Uij=[[2, 0], [2]])
        >>> Rho(H=H, a=1).run()
        KHive(n=3, alpha=[2, 1, 0], beta=[1, 1, 1], gamma=[0, 0, 0], Uij=[[1, 0], [1]])
        """  # noqa: B950
        n: int = self.H.n
        alpha: List[int] = self.compute_alpha()
        beta: List[int] = self.compute_beta()
        gamma: List[int] = [0 for _ in range(self.H.n)]
        Uij: List[List[int]] = self.compute_Uij()

        return KHive(n=n, alpha=alpha, beta=beta, gamma=gamma, Uij=Uij)
