from collections.abc import Callable
from copy import deepcopy
from dataclasses import dataclass

from khive_crystal.khive import KHive
from khive_crystal.khives import KHives


@dataclass(frozen=True)
class FundamentalKHives(KHives):
    """
    This class is a subclass of KHives,
    and has the crystal structure on a set of K-hives with fundamental weight.
    """

    def __post_init__(self) -> None:
        self.is_fundamental_weight()

    def is_fundamental_weight(self) -> None:
        if not all(x in [0, 1] for x in self.alpha):
            raise ValueError("Invaid value! alpha is not a fundamental weight.")

    def phi(self, i: int) -> Callable[[KHive], int]:
        """Compute varphi_i(H)

        Args:
            i (int): i in I.

        Returns:
            Callable: _phi

        Examples:
            >>> H_alpha = FundamentalKHives(n=3, alpha=[1, 1, 0])
            >>> H = H_alpha.highest_weight_vector()
            >>> H_alpha.phi(i=1)(H=H)
            0

            >>> H_alpha = FundamentalKHives(n=3, alpha=[1, 1, 0])
            >>> H = H_alpha.highest_weight_vector()
            >>> H_alpha.phi(i=2)(H=H)
            1
        """

        def _phi(H: KHive) -> int:
            return max(self.inner_product(i=i, weight=self.weight(H=H)), 0)

        return _phi

    def f(self, i: int) -> Callable[[KHive], KHive | None]:
        """Get f_i(H)

        Args:
            i (int): i in I.

        Returns:
            Callable: f_i

        Examples:
            >>> H_alpha = FundamentalKHives(n=3, alpha=[1, 1, 0])
            >>> H = H_alpha.highest_weight_vector()
            >>> H_alpha.f(i=1)(H=H)

            >>> H_alpha = FundamentalKHives(n=3, alpha=[1, 1, 0])
            >>> H = H_alpha.highest_weight_vector()
            >>> H_alpha.f(i=2)(H=H)
            KHive(n=3, alpha=[1, 1, 0], beta=[1, 0, 1], gamma=[0, 0, 0], Uij=[[0, 0], [1]])

            >>> H_alpha = FundamentalKHives(n=3, alpha=[1, 1, 0])
            >>> H: KHive = KHive(n=3, alpha=[1, 1, 0], beta=[1, 0, 1], gamma=[0, 0, 0], Uij=[[0, 0], [1]])
            >>> H_alpha.f(i=1)(H)
            KHive(n=3, alpha=[1, 1, 0], beta=[0, 1, 1], gamma=[0, 0, 0], Uij=[[1, 0], [1]])
        """  # noqa: B950

        def _f(H: KHive) -> KHive | None:
            if self.phi(i=i)(H) == 0:
                return None

            i_as_index: int = i - 1

            beta: list[int] = H.beta.copy()
            gamma: list[int] = [0] * self.n
            Uij: list[list[int]] = deepcopy(H.Uij)

            beta[i_as_index] += -1
            beta[i_as_index + 1] += 1

            act_point_as_index: int
            Uxi: list[int] = H.get_Uji()[i_as_index - 1]
            is_Uxi_geq_zeros: list[bool] = [Uji > 0 for Uji in Uxi]
            if (i == 1) | (not any(is_Uxi_geq_zeros)):
                act_point_as_index = i_as_index
            else:
                act_point_as_index = is_Uxi_geq_zeros.index(True)

            i_as_index_in_Uij: int = i_as_index - act_point_as_index - 1
            ip1_as_index_in_Uij: int = i_as_index - act_point_as_index
            if 0 <= (i_as_index_in_Uij):
                Uij[act_point_as_index][i_as_index_in_Uij] += -1  # U_{ki}
            Uij[act_point_as_index][ip1_as_index_in_Uij] += 1  # U_{k,i+1}

            return KHive(n=self.n, alpha=self.alpha, beta=beta, gamma=gamma, Uij=Uij)

        return _f

    def epsilon(self, i: int) -> Callable[[KHive], int]:
        """Compute varepsilon_i(H)

        Args:
            i (int): i in I.

        Returns:
            Callable: _epsilon

        Examples
            >>> H_alpha = FundamentalKHives(n=3, alpha=[1, 1, 0])
            >>> H = H_alpha.highest_weight_vector()
            >>> H_alpha.epsilon(i=1)(H=H)
            0

            >>> H_alpha = FundamentalKHives(n=3, alpha=[1, 1, 0])
            >>> H = H_alpha.highest_weight_vector()
            >>> H_alpha.epsilon(i=2)(H=H)
            0
        """

        def _epsilon(H: KHive) -> int:
            return max(-self.inner_product(i=i, weight=self.weight(H=H)), 0)

        return _epsilon

    def e(self, i: int) -> Callable[[KHive], KHive | None]:
        """Get e_i(H)

        Args:
            i (int): i in I.

        Returns:
            Callable: e_i

        Examples:
            >>> H_alpha = FundamentalKHives(n=3, alpha=[1, 1, 0])
            >>> H = H_alpha.highest_weight_vector()
            >>> H_alpha.e(i=1)(H=H)

            >>> H_alpha = FundamentalKHives(n=3, alpha=[1, 1, 0])
            >>> H = KHive(
            ...     n=3,
            ...     alpha=[1, 1, 0],
            ...     beta=[1, 0, 1],
            ...     gamma=[0, 0, 0],
            ...     Uij=[[0, 0], [1]]
            ... )
            >>> H_alpha.e(i=2)(H=H)
            KHive(n=3, alpha=[1, 1, 0], beta=[1, 1, 0], gamma=[0, 0, 0], Uij=[[0, 0], [0]])

            >>> H_alpha = FundamentalKHives(n=3, alpha=[1, 1, 0])
            >>> H: KHive = KHive(n=3, alpha=[1, 1, 0], beta=[0, 1, 1], gamma=[0, 0, 0], Uij=[[1, 0], [1]])
            >>> H_alpha.e(i=1)(H)
            KHive(n=3, alpha=[1, 1, 0], beta=[1, 0, 1], gamma=[0, 0, 0], Uij=[[0, 0], [1]])
        """  # noqa: B950

        def _e(H: KHive) -> KHive | None:
            if self.epsilon(i=i)(H) == 0:
                return None

            i_as_index: int = i - 1

            beta: list[int] = H.beta.copy()
            gamma: list[int] = [0] * self.n
            Uij: list[list[int]] = deepcopy(H.Uij)

            beta[i_as_index] += 1
            beta[i_as_index + 1] += -1

            act_point_as_index: int
            Uxip1: list[int] = H.get_Uji()[i_as_index]
            is_Uxip1_geq_zeros: list[bool] = [Uji > 0 for Uji in Uxip1]
            if not any(is_Uxip1_geq_zeros):
                act_point_as_index = i_as_index
            else:
                act_point_as_index = is_Uxip1_geq_zeros.index(True)

            i_as_index_in_Uij: int = i_as_index - act_point_as_index - 1
            ip1_as_index_in_Uij: int = i_as_index - act_point_as_index
            if 0 <= (i_as_index_in_Uij):
                Uij[act_point_as_index][i_as_index_in_Uij] += +1  # U_{ki}
            Uij[act_point_as_index][ip1_as_index_in_Uij] += -1  # U_{k,i+1}

            return KHive(n=self.n, alpha=self.alpha, beta=beta, gamma=gamma, Uij=Uij)

        return _e
