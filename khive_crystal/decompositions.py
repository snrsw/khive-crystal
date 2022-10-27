from copy import deepcopy
from typing import Callable, List, Optional, Union

from khive_crystal.khive import KHive
from khive_crystal.utils import flatten_list


def psi_lambda(H: KHive) -> Union[KHive, List[KHive]]:
    """Split KHive to a pair of KHive and FundamentalKHive.

    Args:
        H (KHive): KHive

    Returns:
        List[KHive]: A pair of KHive and FundamentalKHive.

    Examples:
        >>> H: KHive = KHive(
        ...    n=3,
        ...    alpha=[3, 3, 0],
        ...    beta=[2, 3, 1],
        ...    gamma=[0, 0, 0],
        ...    Uij=[[1, 0], [1]]
        ... )
        >>> psi_lambda(H=H)
        [KHive(n=3, alpha=[2, 2, 0], beta=[1, 2, 1], gamma=[0, 0, 0], Uij=[[1, 0], [1]]), KHive(n=3, alpha=[1, 1, 0], beta=[1, 1, 0], gamma=[0, 0, 0], Uij=[[0, 0], [0]])]
    """  # noqa: B950
    return Split(H=H).run()


class Split:
    """This class has methods to split KHive to a pair of KHive and FundamentalKHive.
    The entry point of this methods is "run". The above function split is a wrapper of run,
    then use split instead of using this class direct.
    """

    def __init__(self, H: KHive) -> None:
        self.H: KHive = H

    def split_alpha(self) -> List[List[int]]:
        """Split H.alpha

        Returns:
            List[List[int]]: [alpha^{(1)}, alpha^{(2)}]

        Examples:
            >>> H: KHive = KHive(
            ...    n=3,
            ...    alpha=[3, 3, 0],
            ...    beta=[2, 3, 1],
            ...    gamma=[0, 0, 0],
            ...    Uij=[[1, 0], [1]]
            ... )
            >>> Split(H=H).split_alpha()
            [[2, 2, 0], [1, 1, 0]]
        """
        alpha: List[int] = self.H.alpha.copy()
        alpha2: List[int] = [1 if alpha_i > 0 else 0 for alpha_i in alpha]
        alpha1: List[int] = [
            alpha_i - alpha2_i for alpha_i, alpha2_i in zip(alpha, alpha2)
        ]

        return [alpha1, alpha2]

    def split_full_Uij(self) -> List[List[List[int]]]:
        """Split H.Uij

        Returns:
            List[List[List[int]]]: [U_{ij}^{(1)}, U_{ij}^{(2)}] (i<=j)

        Examples:
            >>> H: KHive = KHive(
            ...    n=3,
            ...    alpha=[3, 3, 0],
            ...    beta=[2, 3, 1],
            ...    gamma=[0, 0, 0],
            ...    Uij=[[1, 0], [1]]
            ... )
            >>> Split(H=H).split_full_Uij()
            [[[1, 1, 0], [1, 1], [0]], [[1, 0, 0], [1, 0], [0]]]

            >>> H: KHive = KHive(
            ...    n=4,
            ...    alpha=[3, 3, 0, 0],
            ...    beta=[2, 3, 1, 0],
            ...    gamma=[0, 0, 0, 0],
            ...    Uij=[[1, 0, 0], [1, 0], [0]]
            ... )
            >>> Split(H=H).split_full_Uij()
            [[[1, 1, 0, 0], [1, 1, 0], [0, 0], [0]], [[1, 0, 0, 0], [1, 0, 0], [0, 0], [0]]]
        """
        full_Uij: List[List[int]] = deepcopy(self.H.full_Uij)

        # a set of \min\{j in [n] \mid U_{ij} > 0\}s for i in I.
        full_Uij_pos_indexes: List[Optional[int]] = [
            [uij > 0 for uij in Ui].index(True) if sum(Ui) > 0 else None
            for Ui in full_Uij[:-1]
        ]

        full_Uij2: List[List[int]] = [
            [1 if j == Ui_pos_index else 0 for j, uij in enumerate(Ui)]
            for Ui, Ui_pos_index in zip(full_Uij, full_Uij_pos_indexes)
        ] + [
            [0]
        ]  # U_{nn} is always 0
        full_Uij1: List[List[int]] = [
            [uij - uij2 for uij, uij2 in zip(Ui, Ui2)]
            for Ui, Ui2 in zip(full_Uij, full_Uij2)
        ]

        return [full_Uij1, full_Uij2]

    def get_splitted_Uij(self) -> List[List[List[int]]]:
        """Get splitted result of H.Uij from splitted H.full_Uij

        Returns:
            List[List[List[int]]]: [U_{ij}^{(1)}, U_{ij}^{(2)}] (i<j)

        Examples:
            >>> H: KHive = KHive(
            ...    n=3,
            ...    alpha=[3, 3, 0],
            ...    beta=[2, 3, 1],
            ...    gamma=[0, 0, 0],
            ...    Uij=[[1, 0], [1]]
            ... )
            >>> Split(H=H).get_splitted_Uij()
            [[[1, 0], [1]], [[0, 0], [0]]]
        """
        splitted_full_Uij: List[List[List[int]]] = self.split_full_Uij()
        return [
            [full_Ui[1:] for full_Ui in full_Uij[:-1]] for full_Uij in splitted_full_Uij
        ]

    def split_beta(self) -> List[List[int]]:
        """Split H.beta

        Returns:
            List[List[int]]: [alpha^{(1)}, alpha^{(2)}]

        Examples:
            >>> H: KHive = KHive(
            ...    n=3,
            ...    alpha=[3, 3, 0],
            ...    beta=[2, 3, 1],
            ...    gamma=[0, 0, 0],
            ...    Uij=[[1, 0], [1]]
            ... )
            >>> Split(H=H).split_beta()
            [[1, 2, 1], [1, 1, 0]]
        """

        splitted_full_Uij: List[List[List[int]]] = self.split_full_Uij()

        align: Callable[[List[List[int]]], List[List[int]]] = lambda full_Uij: [
            [0] * ((self.H.n) - len(ui)) + ui for ui in full_Uij
        ]
        aligned_full_Uij1: List[List[int]] = align(splitted_full_Uij[0])
        aligned_full_Uij2: List[List[int]] = align(splitted_full_Uij[1])

        get_full_Uji: Callable[
            [List[List[int]]], List[List[int]]
        ] = lambda aligned_full_Uij: [
            list(uj)[: j + 1] for j, uj in enumerate(zip(*aligned_full_Uij))
        ]
        full_Uji1: List[List[int]] = get_full_Uji(aligned_full_Uij1)
        full_Uji2: List[List[int]] = get_full_Uji(aligned_full_Uij2)

        beta1: List[int] = [sum(Uj) for Uj in full_Uji1]
        beta2: List[int] = [sum(Uj) for Uj in full_Uji2]

        return [beta1, beta2]

    def run(self) -> Union[List[KHive], KHive]:
        """Split KHive to a pair of KHive and FundamentalKHive.

        Returns:
            Union[List[KHive], KHive]: A pair of KHive and FundamentalKHive.
        """
        if self.H.is_fundamental_khive():
            return self.H
        splitted_alpha: List[List[int]] = self.split_alpha()
        splitted_Uij: List[List[List[int]]] = self.get_splitted_Uij()
        splitted_beta: List[List[int]] = self.split_beta()
        gamma: List[int] = [0] * self.H.n
        return [
            KHive(
                n=self.H.n,
                alpha=splitted_alpha[i],
                beta=splitted_beta[i],
                gamma=gamma,
                Uij=splitted_Uij[i],
            )
            for i in range(2)
        ]


def psi(H: KHive) -> List[KHive]:
    """Decompose KHive into FundamentalKHive by applying the function split repeatedly.

    Args:
        H (KHive): KHive

    Returns:
        List[KHive]: A pair of KHive and FundamentalKHive.

    Examples:
        >>> H: KHive = KHive(
        ...    n=3,
        ...    alpha=[3, 3, 0],
        ...    beta=[2, 3, 1],
        ...    gamma=[0, 0, 0],
        ...    Uij=[[1, 0], [1]]
        ... )
        >>> psi(H=H)
        [KHive(n=3, alpha=[1, 1, 0], beta=[0, 1, 1], gamma=[0, 0, 0], Uij=[[1, 0], [1]]), KHive(n=3, alpha=[1, 1, 0], beta=[1, 1, 0], gamma=[0, 0, 0], Uij=[[0, 0], [0]]), KHive(n=3, alpha=[1, 1, 0], beta=[1, 1, 0], gamma=[0, 0, 0], Uij=[[0, 0], [0]])]

        >>> H: KHive = KHive(
        ...    n=3,
        ...    alpha=[1, 1, 0],
        ...    beta=[1, 0, 1],
        ...    gamma=[0, 0, 0],
        ...    Uij=[[0, 0], [1]]
        ... )
        >>> psi(H=H)
        [KHive(n=3, alpha=[1, 1, 0], beta=[1, 0, 1], gamma=[0, 0, 0], Uij=[[0, 0], [1]])]
    """  # noqa: B950

    splitted_hive: Union[KHive, List[KHive]] = psi_lambda(H=H)
    if isinstance(splitted_hive, KHive):
        return [splitted_hive]
    return list(flatten_list([psi_lambda(splitted_hive[0]), splitted_hive[1]]))
