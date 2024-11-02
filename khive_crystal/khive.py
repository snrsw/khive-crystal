from dataclasses import dataclass

from khive_crystal.integer_hive_graph import IntegerHiveGraph


@dataclass(frozen=True)
class KHive(IntegerHiveGraph):
    """This class is a subclass of IntegerHiveGraph, and generate a K-hive."""

    def __post_init__(self) -> None:
        super().__post_init__()
        self.is_k_hive()

    def is_eq_edge_sum(self) -> bool:
        """Check sum(alpha) == sum(beta)

        Returns:
            bool: sum(alpha) == sum(beta)

        Examples:
            >>> H: KHive = KHive(
            ...    n=3,
            ...    alpha=[1, 1, 0],
            ...    beta=[1, 1, 0],
            ...    gamma=[0, 0, 0],
            ...    Uij=[[0, 0], [0]]
            ... )
            >>> H.is_eq_edge_sum()
            True
        """
        return sum(self.alpha) == sum(self.beta)

    def is_partition(self) -> bool:
        """Check alpha is a partition.

        Returns:
            bool: If alpha is partition, return True.

        Examples:
            >>> H: KHive = KHive(
            ...    n=3,
            ...    alpha=[1, 1, 0],
            ...    beta=[1, 1, 0],
            ...    gamma=[0, 0, 0],
            ...    Uij=[[0, 0], [0]]
            ... )
            >>> H.is_partition()
            True
        """

        return all(x >= y for x, y in zip(self.alpha, self.alpha[1:], strict=False))

    def is_zero_partition(self) -> bool:
        """Check gamma is 0

        Returns:
            bool: If gamma = (0^n), return True.

        Examples:
            >>> H: KHive = KHive(
            ...    n=3,
            ...    alpha=[1, 1, 0],
            ...    beta=[1, 1, 0],
            ...    gamma=[0, 0, 0],
            ...    Uij=[[0, 0], [0]]
            ... )
            >>> H.is_zero_partition()
            True
        """

        return (sum(self.gamma) == 0) and all(x >= 0 for x in self.gamma)

    def is_Uij_geq_zero(self) -> bool:
        """Check U_{ij} >= 0.

        Returns:
            bool: If U_{ij} >= 0, return True.

        Examples:
            >>> H: KHive = KHive(
            ...    n=3,
            ...    alpha=[1, 1, 0],
            ...    beta=[1, 1, 0],
            ...    gamma=[0, 0, 0],
            ...    Uij=[[0, 0], [0]]
            ... )
            >>> H.is_Uij_geq_zero()
            True
        """
        return all(x >= 0 for x in sum(self.Uij, []))

    def Lij(self, i: int, j: int) -> int:
        """Compute L_{ij} defined by
            L_{ij} = sum_{k=1}^{j-1}U_{ik} - sum_{k=1}^{j}U_{i+1,k}.

        Args:
            i (int): i in [n]
            j (int): j in [i+1, n]_{mathbb{Z}}

        Returns:
            int: L_{ij}

        Raises:
            ValueError: If (i, j) dose not satisfy 1 <= i < j <= n, raise error.

        Examples:
            >>> H: KHive = KHive(
            ...    n=3,
            ...    alpha=[4, 2, 0],
            ...    beta=[2, 1, 3],
            ...    gamma=[0, 0, 0],
            ...    Uij=[[1, 1], [2]]
            ... )
            >>> H.Lij(i=1, j=2)
            2
        """

        if not ((1 <= i <= self.n - 1) or (1 <= j <= self.n)):
            raise ValueError("i and j must be in [n].")

        if i >= j:
            raise ValueError("i and j must be i < j.")

        i_as_index: int = i - 1
        j_as_index: int = j - 1
        Uik: int = sum(self.full_Uij[i_as_index][: j_as_index - i_as_index])
        Uip1k: int = (
            sum(self.full_Uij[i_as_index + 1][: j_as_index - i_as_index - 1])
            if i_as_index + 1 < self.n - 1
            else 0
        )
        return Uik - Uip1k

    def is_Lij_geq_zero(self) -> bool:
        """Check L_{ij} >= 0

        Returns:
            bool: If all L_{ij} >= 0, return True.
        """
        return all(
            self.Lij(i=i + 1, j=j + 1) >= 0
            for i in range(self.n - 1)
            for j in range(i + 1, self.n)
        )

    def is_Uii_geq_zero(self) -> bool:
        return all(self.compute_Uii(i=i + 1) >= 0 for i in range(self.n))

    def is_k_hive(self) -> None:
        """Check given values saitisfy the definition of K-hive.

        Raises:
            ValueError: _description_
        """
        is_k_hive: bool = all(
            [
                self.is_eq_edge_sum(),
                self.is_partition(),
                self.is_Uij_geq_zero(),
                self.is_Lij_geq_zero(),
                self.is_Uii_geq_zero(),
            ]
        )
        if not is_k_hive:
            raise ValueError("Given values dose not satisfy the definition of K-hive")

    def is_fundamental_khive(self) -> bool:
        """Check self is a fundamental khive.

        Returns:
            bool: If alpha is a fundamental khive, return True
        """
        return all(alpha_i in [0, 1] for alpha_i in self.alpha)


def khive(
    n: int, alpha: list[int], beta: list[int], gamma: list[int], Uij: list[list[int]]
) -> KHive:
    """Get a instance of KHive

    Args:
        n (int): size
        alpha (List[int]): right edge labels
        beta (List[int]): below edge labels
        gamma (List[int]): left edge labels
        Uij (List[List[int]]): rhombi

    Returns:
        KHive: K-hive
    """
    return KHive(n=n, alpha=alpha, beta=beta, gamma=gamma, Uij=Uij)
