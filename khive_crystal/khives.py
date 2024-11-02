from collections.abc import Callable
from copy import deepcopy
from dataclasses import dataclass

from khive_crystal.khive import KHive


@dataclass(frozen=True)
class KHives:
    """This class has the crystal structure on K-hives."""

    n: int
    alpha: list[int]

    def __post_init__(self) -> None:
        pass

    def isin(self, H: KHive) -> bool:
        """Check a given KHive is in a instance.

        Args:
            H (KHive): KHive

        Returns:
            bool: If KHive in KHives, return True, otherwise return False.

        Examples:
            >>> H: KHive = KHive(n=3, alpha=[2, 1, 0], beta=[2, 1, 0], gamma=[0, 0, 0], Uij=[[0, 0], [0]])
            >>> KHives(n=3, alpha=[2, 1, 0]).isin(H)
            True
        """  # noqa: B950
        return (H.alpha == self.alpha) and (H.n == self.n)

    def highest_weight_vector(self) -> KHive:
        """Get the highest weight vector H_{alpha} in mathbb{H}(alpha).

        Returns:
            KHive: H_{alpha} = (alpha, alpha, 0, (0))

        Examples:
            >>> H_alpha = KHives(n=3, alpha=[3, 2, 0])
            >>> H_alpha.highest_weight_vector()
            KHive(n=3, alpha=[3, 2, 0], beta=[3, 2, 0], gamma=[0, 0, 0], Uij=[[0, 0], [0]])
        """
        zeros: list[int] = [0] * self.n
        Uij: list[list[int]] = [[0] * (self.n - i) for i in range(1, self.n)]
        return KHive(n=self.n, alpha=self.alpha, beta=self.alpha, gamma=zeros, Uij=Uij)

    def weight(self, H: KHive) -> list[int]:
        """Get wt(H)

        Args:
            H (KHive): H in mathbb{H}(alpha)

        Returns:
            List[int]: wt(H).

        Examples:
            >>> H_alpha = KHives(n=3, alpha=[3, 2, 0])
            >>> H = H_alpha.highest_weight_vector()
            >>> H_alpha.weight(H=H)
            [3, 2, 0]
        """
        return H.beta

    def inner_product(self, i: int, weight: list[int]) -> int:
        """Compute <h_i, wt(H)>

        Args:
            i (int): h_i
            weight (List[int]): weight.

        Returns:
            int: <h_i, wt(H)>

        Examples:
            >>> H_alpha = KHives(n=3, alpha=[3, 2, 0])
            >>> H = H_alpha.highest_weight_vector()
            >>> H_alpha.inner_product(i=1, weight=H_alpha.weight(H=H))
            1
        """
        if i not in [_ + 1 for _ in range(self.n - 1)]:
            raise ValueError("i must be in I.")

        i_as_index: int = i - 1
        return weight[i_as_index] - weight[i_as_index + 1]

    def f(self, i: int) -> Callable[[KHive], KHive | None]:
        """Get f_i(H)

        Args:
            i (int): i in I,

        Returns:
            Callable: f_i

        Examples:
            >>> H_alpha = KHives(n=3, alpha=[2, 2, 0])
            >>> H = H_alpha.highest_weight_vector()
            >>> H_alpha.f(i=1)(H=H)

            >>> H_alpha = KHives(n=3, alpha=[3, 2, 0])
            >>> H = H_alpha.highest_weight_vector()
            >>> H_alpha.f(i=1)(H=H)
            KHive(n=3, alpha=[3, 2, 0], beta=[2, 3, 0], gamma=[0, 0, 0], Uij=[[1, 0], [0]])
        """

        def f_i(H: KHive) -> KHive | None:
            if self.phi(i=i)(H) == 0:
                return None

            beta: list[int] = deepcopy(H.beta)
            Uij: list[list[int]] = deepcopy(H.Uij)

            i_as_index: int = i - 1

            beta[i_as_index] += -1
            beta[i_as_index + 1] += 1

            act_point_dual: int = [
                self._phi(i=i, k=k)(H) > 0 for k in range(i, -1, -1)
            ].index(False) - 1
            act_point: int = i - act_point_dual
            act_point_as_index: int = act_point - 1

            i_as_index_in_Uij: int = i_as_index - act_point_as_index - 1
            ip1_as_index_in_Uij: int = i_as_index - act_point_as_index
            if 0 <= (i_as_index_in_Uij):
                Uij[act_point_as_index][i_as_index_in_Uij] += -1  # U_{ki}
            Uij[act_point_as_index][ip1_as_index_in_Uij] += 1  # U_{k,i+1}

            return KHive(n=H.n, alpha=H.alpha, beta=beta, gamma=H.gamma, Uij=Uij)

        return f_i

    def phi(self, i: int) -> Callable[[KHive], int]:
        """Compute varphi_i(H)

        Args:
            i (int): i in I.

        Returns:
            Callable: varphi_i

        Examples:
            >>> H_alpha = KHives(n=3, alpha=[3, 2, 0])
            >>> H = H_alpha.highest_weight_vector()
            >>> H_alpha.phi(i=1)(H=H)
            1
        """

        def _phi(H: KHive) -> int:
            return self._phi(i=i, k=i)(H)

        return _phi

    def _phi(self, i: int, k: int) -> Callable[[KHive], int]:
        """Compute varphi_i^k(H)

        Args:
            i (int): i in I.
            k (int): k in [i].

        Returns:
            Callable: varphi_i^j

        Examples:
            >>> H_alpha = KHives(n=3, alpha=[3, 2, 0])
            >>> H = H_alpha.highest_weight_vector()
            >>> H_alpha._phi(i=1, k=1)(H=H)
            1

            >>> H_alpha = KHives(n=3, alpha=[3, 2, 0])
            >>> H = H_alpha.highest_weight_vector()
            >>> H_alpha._phi(i=2, k=1)(H=H)
            0

            >>> H_alpha = KHives(n=3, alpha=[3, 2, 0])
            >>> H = H_alpha.highest_weight_vector()
            >>> H_alpha._phi(i=2, k=2)(H=H)
            2
        """

        def phi_i_k(H: KHive) -> int:
            i_as_index: int = i - 1
            k_as_index: int = k - 1

            if k == 0:
                return 0

            hat_phi_i_k = (
                self._phi(i=i, k=k - 1)(H)
                + H.full_Uij[k_as_index][i_as_index - k_as_index]  # U_{ki}
                - H.full_Uij[k_as_index + 1][i_as_index - k_as_index]  # U_{k+1,i+1}
            )

            return max(hat_phi_i_k, 0)

        return phi_i_k

    def e(self, i: int) -> Callable[[KHive], KHive | None]:
        """Get e_i(H)

        Args:
            i (int): i in I,

        Returns:
            Callable: e_i

        Examples:
            >>> H_alpha = KHives(n=3, alpha=[2, 2, 0])
            >>> H = H_alpha.highest_weight_vector()
            >>> H_alpha.e(i=1)(H=H)

            >>> H_alpha = KHives(n=3, alpha=[3, 2, 0])
            >>> H: KHive = KHive(
            ...    n=3,
            ...    alpha=[3, 2, 0],
            ...    beta=[2, 3, 0],
            ...    gamma=[0, 0, 0],
            ...    Uij=[[1, 0], [0]]
            ... )
            >>> H_alpha.e(i=1)(H=H)
            KHive(n=3, alpha=[3, 2, 0], beta=[3, 2, 0], gamma=[0, 0, 0], Uij=[[0, 0], [0]])
        """

        def e_i(H: KHive) -> KHive | None:
            if self.epsilon(i=i)(H) == 0:
                return None

            beta: list[int] = deepcopy(H.beta)
            Uij: list[list[int]] = deepcopy(H.Uij)

            i_as_index: int = i - 1

            beta[i_as_index] += 1
            beta[i_as_index + 1] += -1

            act_point_dual: int = [
                self._epsilon(i=i, k=k)(H) > 0 for k in range(i + 1, -1, -1)
            ].index(False) - 1
            act_point: int = i + 1 - (i - act_point_dual)
            act_point_as_index: int = act_point - 1

            i_as_index_in_Uij: int = i_as_index - act_point_as_index - 1
            ip1_as_index_in_Uij: int = i_as_index - act_point_as_index
            if 0 <= (i_as_index_in_Uij):
                Uij[act_point_as_index][i_as_index_in_Uij] += +1  # U_{ki}
            Uij[act_point_as_index][ip1_as_index_in_Uij] += -1  # U_{k,i+1}

            return KHive(n=H.n, alpha=H.alpha, beta=beta, gamma=H.gamma, Uij=Uij)

        return e_i

    def epsilon(self, i: int) -> Callable[[KHive], int]:
        """Compute vare[silon_i(H)

        Args:
            i (int): i in I.

        Returns:
            Callable: varepsilon_i

        Examples:
            >>> H_alpha = KHives(n=3, alpha=[3, 2, 0])
            >>> H = H_alpha.highest_weight_vector()
            >>> H_alpha.epsilon(i=1)(H=H)
            0

            >>> H_alpha = KHives(n=3, alpha=[3, 2, 0])
            >>> H: KHive = KHive(
            ...    n=3,
            ...    alpha=[3, 2, 0],
            ...    beta=[2, 3, 0],
            ...    gamma=[0, 0, 0],
            ...    Uij=[[1, 0], [0]]
            ... )
            >>> H_alpha.epsilon(i=1)(H=H)
            1
        """

        def _epsilon(H: KHive) -> int:
            return self._epsilon(i=i, k=i + 1)(H)

        return _epsilon

    def _epsilon(self, i: int, k: int) -> Callable[[KHive], int]:
        """Compute varepsilon_i^k(H)

        Args:
            i (int): i in I.
            k (int): k in [i+1].

        Returns:
            Callable: varepsilon_i^j

        Examples:
            >>> H_alpha = KHives(n=3, alpha=[3, 2, 0])
            >>> H = H_alpha.highest_weight_vector()
            >>> H_alpha._epsilon(i=1, k=1)(H=H)
            0

            >>> H_alpha = KHives(n=3, alpha=[3, 2, 0])
            >>> H = H_alpha.highest_weight_vector()
            >>> H_alpha._epsilon(i=2, k=1)(H=H)
            0

            >>> H_alpha = KHives(n=3, alpha=[3, 2, 0])
            >>> H: KHive = KHive(
            ...    n=3,
            ...    alpha=[3, 2, 0],
            ...    beta=[2, 3, 0],
            ...    gamma=[0, 0, 0],
            ...    Uij=[[1, 0], [0]]
            ... )
            >>> H_alpha._epsilon(i=1, k=1)(H=H)
            0

            >>> H_alpha = KHives(n=3, alpha=[3, 2, 0])
            >>> H: KHive = KHive(
            ...     n=3,
            ...     alpha=[3, 2, 0],
            ...     beta=[2, 3, 0],
            ...     gamma=[0, 0, 0],
            ...     Uij=[[1, 0], [0]]
            ... )
            >>> H_alpha._epsilon(i=1, k=2)(H=H)
            1
        """

        def epsilon_i_k(H: KHive) -> int:
            i_as_index: int = i - 1

            if k == 0:
                return 0

            hat_epsilon_i_k: int
            if k == i + 1:
                hat_epsilon_i_k = (
                    self._epsilon(i=i, k=k - 1)(H)
                    + H.full_Uij[i_as_index + 2 - k][
                        i_as_index - (i_as_index + 1 - k)
                    ]  # U_{i+2-k,i+1}
                )
            else:
                hat_epsilon_i_k = (
                    self._epsilon(i=i, k=k - 1)(H)
                    + H.full_Uij[i_as_index + 2 - k][
                        i_as_index - (i_as_index + 1 - k)
                    ]  # U_{i+2-k,i+1}
                    - H.full_Uij[i_as_index + 1 - k][
                        i_as_index - (i_as_index + 1 - k)
                    ]  # U_{i+1-k,i}
                )

            return max(hat_epsilon_i_k, 0)

        return epsilon_i_k


def khives(n: int, alpha: list[int]) -> KHives:
    """Get a instance of KHives

    Args:
        n (int): size
        alpha (List[int]): right edge labels

    Returns:
        KHives: Khives
    """
    return KHives(n=n, alpha=alpha)
