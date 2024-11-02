from collections.abc import Callable
from copy import deepcopy

from khive_crystal.khive import KHive
from khive_crystal.utils import flatten_list


def psi_lambda(H: KHive) -> KHive | list[KHive]:
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

    def split_alpha(self) -> list[list[int]]:
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
        alpha: list[int] = self.H.alpha.copy()
        alpha2: list[int] = [1 if alpha_i > 0 else 0 for alpha_i in alpha]
        alpha1: list[int] = [
            alpha_i - alpha2_i for alpha_i, alpha2_i in zip(alpha, alpha2, strict=False)
        ]

        return [alpha1, alpha2]

    def split_full_Uij(self) -> list[list[list[int]]]:
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
        full_Uij: list[list[int]] = deepcopy(self.H.full_Uij)

        # a set of \min\{j in [n] \mid U_{ij} > 0\}s for i in I.
        full_Uij_pos_indexes: list[int | None] = [
            [uij > 0 for uij in Ui].index(True) if sum(Ui) > 0 else None
            for Ui in full_Uij[:-1]
        ]

        full_Uij2: list[list[int]] = [
            [1 if j == Ui_pos_index else 0 for j, uij in enumerate(Ui)]
            for Ui, Ui_pos_index in zip(full_Uij, full_Uij_pos_indexes, strict=False)
        ] + [
            [0]
        ]  # U_{nn} is always 0
        full_Uij1: list[list[int]] = [
            [uij - uij2 for uij, uij2 in zip(Ui, Ui2, strict=False)]
            for Ui, Ui2 in zip(full_Uij, full_Uij2, strict=False)
        ]

        return [full_Uij1, full_Uij2]

    def get_splitted_Uij(self) -> list[list[list[int]]]:
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
        splitted_full_Uij: list[list[list[int]]] = self.split_full_Uij()
        return [
            [full_Ui[1:] for full_Ui in full_Uij[:-1]] for full_Uij in splitted_full_Uij
        ]

    def split_beta(self) -> list[list[int]]:
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

        splitted_full_Uij: list[list[list[int]]] = self.split_full_Uij()

        align: Callable[[list[list[int]]], list[list[int]]] = lambda full_Uij: [
            [0] * ((self.H.n) - len(ui)) + ui for ui in full_Uij
        ]
        aligned_full_Uij1: list[list[int]] = align(splitted_full_Uij[0])
        aligned_full_Uij2: list[list[int]] = align(splitted_full_Uij[1])

        get_full_Uji: Callable[
            [list[list[int]]], list[list[int]]
        ] = lambda aligned_full_Uij: [
            list(uj)[: j + 1] for j, uj in enumerate(zip(*aligned_full_Uij, strict=False))
        ]
        full_Uji1: list[list[int]] = get_full_Uji(aligned_full_Uij1)
        full_Uji2: list[list[int]] = get_full_Uji(aligned_full_Uij2)

        beta1: list[int] = [sum(Uj) for Uj in full_Uji1]
        beta2: list[int] = [sum(Uj) for Uj in full_Uji2]

        return [beta1, beta2]

    def run(self) -> list[KHive] | KHive:
        """Split KHive to a pair of KHive and FundamentalKHive.

        Returns:
            Union[List[KHive], KHive]: A pair of KHive and FundamentalKHive.
        """
        if self.H.is_fundamental_khive():
            return self.H
        splitted_alpha: list[list[int]] = self.split_alpha()
        splitted_Uij: list[list[list[int]]] = self.get_splitted_Uij()
        splitted_beta: list[list[int]] = self.split_beta()
        gamma: list[int] = [0] * self.H.n
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


def psi(H: KHive) -> list[KHive]:
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
        ...    alpha=[3, 1, 0],
        ...    beta=[2, 2, 0],
        ...    gamma=[0, 0, 0],
        ...    Uij=[[1, 0], [0]]
        ... )
        >>> psi(H=H)
        [KHive(n=3, alpha=[1, 0, 0], beta=[0, 1, 0], gamma=[0, 0, 0], Uij=[[1, 0], [0]]), KHive(n=3, alpha=[1, 0, 0], beta=[1, 0, 0], gamma=[0, 0, 0], Uij=[[0, 0], [0]]), KHive(n=3, alpha=[1, 1, 0], beta=[1, 1, 0], gamma=[0, 0, 0], Uij=[[0, 0], [0]])]

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

    splitted_hive: KHive | list[KHive] = psi_lambda(H=H)
    if isinstance(splitted_hive, KHive):
        return [splitted_hive]
    return list(flatten_list([psi(splitted_hive[0]), splitted_hive[1]]))


def psi_inv(H: list[KHive]) -> KHive:
    """Compose tensor product of K-hives to a K-hive by Psi^{-1}.

    Args:
        H (List[KHive]): list of K-hives expected to belong to an image of Psi.

    Returns:
        KHive: Psi^{-1}(H)

    Examples:
        >>> H1: KHive = KHive(n=3, alpha=[1, 1, 0], beta=[1, 1, 0], gamma=[0, 0, 0], Uij=[[0, 0], [0]])
        >>> H2: KHive = KHive(n=3, alpha=[1, 1, 0], beta=[0, 1, 1], gamma=[0, 0, 0], Uij=[[1, 0], [1]])
        >>> H: List[KHive] = [H1, H2]
        >>> psi_inv(H=H)
        KHive(n=3, alpha=[2, 2, 0], beta=[1, 2, 1], gamma=[0, 0, 0], Uij=[[1, 0], [1]])
    """  # noqa: B950
    return Compose(H=H).run()


class Compose:
    """This class has methods to compose tensor products of KHive to a KHive.
    The entry point of this methods is "run". The function psi_inv is a wrapper of run,
    then use psi_inv instead of using this class direct.
    """

    def __init__(self, H: list[KHive]) -> None:
        self.H: list[KHive] = H
        self.validate_is_khive()
        self.validate_khive_size()

    def validate_is_khive(self) -> None:
        """Validate that H is a list of KHive.

        Raises:
            ValueError: If H has an object which is not a Khive, raise error.

        Examples:
            >>> H1: KHive = KHive(n=3, alpha=[1, 1, 0], beta=[1, 1, 0], gamma=[0, 0, 0], Uij=[[0, 0], [0]])
            >>> H: List[int] = [H1, 1]
            >>> Compose(H=H)
            Traceback (most recent call last):
            ...
            ValueError: H must be a list of KHive!
        """  # noqa: B950
        is_valid: bool = all(isinstance(Hi, KHive) for Hi in self.H)
        if not is_valid:
            raise ValueError("H must be a list of KHive!")

    def validate_khive_size(self) -> None:
        """Validate that H.n

        Raises:
            ValueError: If H has an KHive with n != H[0].n. raise error.

        Examples:
            >>> H1: KHive = KHive(n=3, alpha=[1, 1, 0], beta=[1, 1, 0], gamma=[0, 0, 0], Uij=[[0, 0], [0]])
            >>> H2: KHive = KHive(n=2, alpha=[1, 0], beta=[1, 0], gamma=[0, 0], Uij=[[0]])
            >>> H: List[int] = [H1, H2]
            >>> Compose(H=H)
            Traceback (most recent call last):
            ...
            ValueError: All elements of H must have the same n
        """  # noqa: B950
        lead_n: int = self.H[0].n
        is_valid: bool = all(Hi.n == lead_n for Hi in self.H)
        if not is_valid:
            raise ValueError("All elements of H must have the same n")

    def get_n(self) -> int:
        """Get n of H

        Returns:
            int: H[0].n
        """
        return self.H[0].n

    def compose_alpha(self) -> list[int]:
        """Compose alpha_i

        Returns:
            List[int]: alpha of psi^{-1}(H)

        Examples:
            >>> H1: KHive = KHive(n=3, alpha=[1, 1, 0], beta=[1, 1, 0], gamma=[0, 0, 0], Uij=[[0, 0], [0]])
            >>> H: List[int] = [H1, H1]
            >>> Compose(H=H).compose_alpha()
            [2, 2, 0]
        """  # noqa: B950
        return [sum(alpha_i) for alpha_i in zip(*[Hi.alpha for Hi in self.H], strict=False)]

    def compose_beta(self) -> list[int]:
        """Compose beta_i

        Returns:
            List[int]: beta of psi^{-1}(H)

        Examples:
            >>> H1: KHive = KHive(n=3, alpha=[1, 1, 0], beta=[1, 1, 0], gamma=[0, 0, 0], Uij=[[0, 0], [0]])
            >>> H2: KHive = KHive(n=3, alpha=[1, 1, 0], beta=[0, 1, 1], gamma=[0, 0, 0], Uij=[[1, 0], [1]])
            >>> H: List[int] = [H1, H2]
            >>> Compose(H=H).compose_beta()
            [1, 2, 1]
        """  # noqa: B950
        return [sum(beta_i) for beta_i in zip(*[Hi.beta for Hi in self.H], strict=False)]

    def compose_Uij(self) -> list[list[int]]:
        """Compose U_{ij}^{k}

        Returns:
            List[int]: U_{ij} of psi^{-1}(H)

        Examples:
            >>> H1: KHive = KHive(n=3, alpha=[1, 1, 0], beta=[1, 1, 0], gamma=[0, 0, 0], Uij=[[0, 0], [0]])
            >>> H: List[int] = [H1, H1]
            >>> Compose(H=H).compose_Uij()
            [[0, 0], [0]]

            >>> H1: KHive = KHive(n=3, alpha=[1, 1, 0], beta=[1, 1, 0], gamma=[0, 0, 0], Uij=[[0, 0], [0]])
            >>> H2: KHive = KHive(n=3, alpha=[1, 1, 0], beta=[0, 1, 1], gamma=[0, 0, 0], Uij=[[1, 0], [1]])
            >>> H: List[int] = [H1, H2]
            >>> Compose(H=H).compose_Uij()
            [[1, 0], [1]]

            >>> H2: KHive = KHive(n=3, alpha=[1, 1, 0], beta=[0, 1, 1], gamma=[0, 0, 0], Uij=[[1, 0], [1]])
            >>> H: List[int] = [H2, H2]
            >>> Compose(H=H).compose_Uij()
            [[2, 0], [2]]
        """  # noqa: B950
        return [
            [sum(uij) for uij in zip(*Ui, strict=False)] for Ui in zip(*[Hi.Uij for Hi in self.H], strict=False)
        ]

    def run(self) -> KHive:
        """Generate psi^{-1}

        Returns:
            KHive: psi^{-1}(H)

        Examples:
            >>> H1: KHive = KHive(n=3, alpha=[1, 1, 0], beta=[1, 1, 0], gamma=[0, 0, 0], Uij=[[0, 0], [0]])
            >>> H2: KHive = KHive(n=3, alpha=[1, 1, 0], beta=[0, 1, 1], gamma=[0, 0, 0], Uij=[[1, 0], [1]])
            >>> H: List[KHive] = [H1, H2]
            >>> Compose(H=H).run()
            KHive(n=3, alpha=[2, 2, 0], beta=[1, 2, 1], gamma=[0, 0, 0], Uij=[[1, 0], [1]])
        """  # noqa: B950
        n: int = self.get_n()
        alpha: list[int] = self.compose_alpha()
        beta: list[int] = self.compose_beta()
        gamma: list[int] = [0 for _ in range(n)]
        Uij: list[list[int]] = self.compose_Uij()
        return KHive(n=n, alpha=alpha, beta=beta, gamma=gamma, Uij=Uij)
