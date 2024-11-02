from dataclasses import dataclass, field


def align_length(n: int, edge: list[int]) -> list[int]:
    """align length of edge to n

    Args:
        n (int): length
        edge (List[int]): boundary edge labels

    Returns:
        List[int]: length aligned boundary edge labels

    Examples:
        >>> align_length(n=3, edge=[2, 2])
        [2, 2, 0]
    """
    return edge if len(edge) == n else edge + [0] * (n - len(edge))


@dataclass(frozen=True)
class IntegerHiveGraph:
    """
    This class generate a integer hive graph.
    """

    n: int
    alpha: list[int]
    beta: list[int]
    gamma: list[int]
    Uij: list[list[int]]
    full_Uij: list[list[int]] = field(init=False, repr=False)

    def __post_init__(self) -> None:
        object.__setattr__(self, "alpha", align_length(n=self.n, edge=self.alpha))
        object.__setattr__(self, "beta", align_length(n=self.n, edge=self.beta))
        object.__setattr__(self, "gamma", align_length(n=self.n, edge=self.gamma))
        object.__setattr__(self, "full_Uij", self.add_Uii())
        self.is_integer_hive_graph()

    def compute_Uii(self, i: int) -> int:
        """compute U_{ii}
            U_{ii} := beta_{i} - sum_{k=1}^{i-1}U_{ki}

        Args:
            i (int): index

        Returns:
            int: U_{ii}

        Examples:
            >>> H: IntegerHiveGraph = IntegerHiveGraph(
            ...    n=3,
            ...    alpha=[1, 1, 0],
            ...    beta=[1, 1, 0],
            ...    gamma=[0, 0, 0],
            ...    Uij=[[0, 0], [0]]
            ... )
            >>> H.compute_Uii(i=1)
            1

            >>> H: IntegerHiveGraph = IntegerHiveGraph(
            ...    n=3,
            ...    alpha=[3, 3, 0],
            ...    beta=[2, 3, 1],
            ...    gamma=[0, 0, 0],
            ...    Uij=[[1, 0], [1]]
            ... )
            >>> H.compute_Uii(i=2)
            2

            >>> H: IntegerHiveGraph = IntegerHiveGraph(
            ...    n=3,
            ...    alpha=[3, 1, 0],
            ...    beta=[1, 1, 2],
            ...    gamma=[0, 0, 0],
            ...    Uij=[[1, 1], [1]]
            ... )
            >>> H.compute_Uii(i=3)
            0
        """

        if i < 1:
            raise ValueError("i should be in [n].")

        i_as_index: int = i - 1
        Uji: list[list[int]] = self.get_Uji()
        return (
            self.beta[i_as_index] - sum(Uji[i_as_index - 1]) if 1 < i else self.beta[0]
        )

    def add_Uii(self) -> list[list[int]]:
        """Add Uii to Uij

        Returns:
            List[List[int]]: [U_{ii}, U_{i,i+1}, ...]

        Examples:
            >>> H: IntegerHiveGraph = IntegerHiveGraph(
            ...    n=3,
            ...    alpha=[1, 1, 0],
            ...    beta=[1, 1, 0],
            ...    gamma=[0, 0, 0],
            ...    Uij=[[0, 0], [0]]
            ... )
            >>> H.add_Uii()
            [[1, 0, 0], [1, 0], [0]]

            >>> H: IntegerHiveGraph = IntegerHiveGraph(
            ...    n=3,
            ...    alpha=[3, 2, 0],
            ...    beta=[3, 2, 0],
            ...    gamma=[0, 0, 0],
            ...    Uij=[[0, 0], [0]]
            ... )
            >>> H.add_Uii()
            [[3, 0, 0], [2, 0], [0]]
        """
        return [[self.compute_Uii(i=i + 1)] + ui for i, ui in enumerate(self.Uij)] + [
            [self.compute_Uii(i=self.n)]
        ]

    def get_Uji(self) -> list[list[int]]:
        """get Uji by aligning Uij.
        Let Uij = [[U_{12}, U_{13}, U_{14}], [U_{23}, U_{24}], [U_{34}]].
        Then Uji = [[U_{12}], [U_{13}, U_{23}], [U_{14}, U_{24}], [U_{34}]].

        Returns:
            List[List[int]]: Uji

        Examples:
            >>> H: IntegerHiveGraph = IntegerHiveGraph(
            ...    n=3,
            ...    alpha=[2, 1, 0],
            ...    beta=[2, 1, 0],
            ...    gamma=[0, 0, 0],
            ...    Uij=[[0, 0], [0]]
            ... )
            >>> H.get_Uji()
            [[0], [0, 0]]

            >>> H: IntegerHiveGraph = IntegerHiveGraph(
            ...    n=3,
            ...    alpha=[2, 1, 0],
            ...    beta=[0, 2, 1],
            ...    gamma=[0, 0, 0],
            ...    Uij=[[2, 0], [1]]
            ... )
            >>> H.get_Uji()
            [[2], [0, 1]]

            >>> H: IntegerHiveGraph = IntegerHiveGraph(
            ...    n=3,
            ...    alpha=[3, 1, 0],
            ...    beta=[1, 1, 2],
            ...    gamma=[0, 0, 0],
            ...    Uij=[[1, 1], [1]]
            ... )
            >>> H.get_Uji()
            [[1], [1, 1]]

            >>> H: IntegerHiveGraph = IntegerHiveGraph(
            ...    n=3,
            ...    alpha=[1, 1, 0],
            ...    beta=[1, 0, 1],
            ...    gamma=[0, 0, 0],
            ...    Uij=[[0, 0], [1]]
            ... )
            >>> H.get_Uji()
            [[0], [0, 1]]
        """

        aligned_Uij: list[list[int]] = [
            [0] * ((self.n - 1) - len(ui)) + ui for ui in self.Uij
        ]
        return [list(uj)[: j + 1] for j, uj in enumerate(zip(*aligned_Uij, strict=False))]

    def integer_hive_graph_condition(self, k: int) -> bool:
        """Compute the condition of integer hive graph:
        beta_{i} = gamma_{i} + sum_{k=1}^{i-1} U_{ki} + alpha_{i} - sum_{k=i+1}^{n}U_{ik}

        Args:
            k (int): i

        Returns:
            bool: is condition holds
        """

        k_as_index: int = k - 1
        is_hold: bool
        if k == 1:
            is_hold = self.beta[0] == (self.gamma[0] + self.alpha[0] - sum(self.Uij[0]))
        elif k == self.n:
            is_hold = self.beta[k_as_index] == (
                self.gamma[k_as_index]
                + sum(self.get_Uji()[k_as_index - 1])  # \sum_{l=1}^{k-1} U_{lk}
                + self.alpha[k_as_index]
            )
        else:
            is_hold = self.beta[k_as_index] == (
                self.gamma[k_as_index]
                + sum(self.get_Uji()[k_as_index - 1])  # \sum_{l=1}^{k-1} U_{lk}
                + self.alpha[k_as_index]
                - sum(self.Uij[k_as_index])
            )  # \sum_{l=k+1}^{n} U_{kl}
        return is_hold

    def is_integer_hive_graph(self) -> None:
        """Check the condition of integer hive graph.

        Raises:
            ValueError: If the condition is not satisfied raise error, otherwise nothing.
        """

        conditions: list[bool] = [
            self.integer_hive_graph_condition(k=k + 1) for k in range(self.n)
        ]

        is_condition_hold: bool = sum(conditions) == len(conditions)

        if not is_condition_hold:
            raise ValueError(
                f"""
                Given values dose not satisfy the definition of a interger hive graph.
                {self},
                conditions = {conditions}
                """
            )
