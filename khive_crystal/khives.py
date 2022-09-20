from copy import deepcopy
from dataclasses import dataclass
from logging import DEBUG, StreamHandler, getLogger
from typing import Callable, List, Optional, Union

from khive_crystal.khive import KHive

logger = getLogger(__name__)
handler = StreamHandler()
handler.setLevel(DEBUG)
logger.setLevel(DEBUG)
logger.addHandler(handler)
logger.propagate = False


@dataclass(frozen=True)
class KHives:
    n: int
    alpha: List[int]

    def __post_init__(self) -> None:
        pass

    def highest_weight_vector(self) -> KHive:
        """Get the highest weight vector H_{alpha} in mathbb{H}(alpha).

        Returns:
            KHive: H_{alpha} = (alpha, alpha, 0, (0))

        Examples:
            >>> H_alpha = KHives(n=3, alpha=[3, 2, 0])
            >>> H_alpha.highest_weight_vector()
            KHive(n=3, alpha=[3, 2, 0], beta=[3, 2, 0], gamma=[0, 0, 0], Uij=[[0, 0], [0]])
        """
        zeros: List[int] = [0] * self.n
        Uij: List[List[int]] = [[0] * (self.n - i) for i in range(1, self.n)]
        return KHive(n=self.n, alpha=self.alpha, beta=self.alpha, gamma=zeros, Uij=Uij)

    def weight(self, H: KHive) -> Callable:
        """Get wt(H) or compute <h_i, wt(H)>

        Args:
            H (KHive): H in mathbb{H}(alpha)

        Returns:
            Callable: If i is given, return <h_i, wt(H)>, otherwise return wt(H).

        Examples:
            >>> H_alpha = KHives(n=3, alpha=[3, 2, 0])
            >>> H = H_alpha.highest_weight_vector()
            >>> H_alpha.weight(H=H)()
            [3, 2, 0]

            >>> H_alpha = KHives(n=3, alpha=[3, 2, 0])
            >>> H = H_alpha.highest_weight_vector()
            >>> H_alpha.weight(H=H)(i=1)
            1
        """

        def _weight(i: Optional[int] = None) -> Union[List[int], int]:
            weight: List[int] = H.beta
            result: Union[List[int], int]
            if i is None:
                result = weight
            else:
                if not (i in [_ + 1 for _ in range(self.n - 1)]):
                    raise ValueError("Invalid value! i must be in I.")

                i_as_index: int = i - 1
                result = weight[i_as_index] - weight[i_as_index + 1]
            return result

        return _weight

    def f(self, i: int) -> Callable:
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

        def f_i(H: KHive) -> Union[KHive, None]:
            if self.phi(i=i)(H=H) == 0:
                return None

            beta: List[int] = deepcopy(H.beta)
            Uij: List[List[int]] = deepcopy(H.Uij)

            i_as_index: int = i - 1

            beta[i_as_index] += -1
            beta[i_as_index + 1] += 1

            act_point_dual: int = [
                self._phi(i=i, k=k)(H=H) > 0 for k in range(i, -1, -1)
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

    def phi(self, i: int) -> Callable:
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
            return self._phi(i=i, k=i)(H=H)

        return _phi

    def _phi(self, i: int, k: int) -> Callable:
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

    def e(self, i: int):
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
            >>> H: KHive = KHive(n=3, alpha=[3, 2, 0], beta=[2, 3, 0], gamma=[0, 0, 0], Uij=[[1, 0], [0]])
            >>> H_alpha.e(i=1)(H=H)
            KHive(n=3, alpha=[3, 2, 0], beta=[3, 2, 0], gamma=[0, 0, 0], Uij=[[0, 0], [0]])
        """

        def e_i(H: KHive) -> Union[KHive, None]:
            if self.epsilon(i=i)(H=H) == 0:
                return None

            beta: List[int] = deepcopy(H.beta)
            Uij: List[List[int]] = deepcopy(H.Uij)

            i_as_index: int = i - 1

            beta[i_as_index] += 1
            beta[i_as_index + 1] += -1

            act_point_dual: int = [
                self._epsilon(i=i, k=k)(H=H) > 0 for k in range(i + 1, -1, -1)
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

    def epsilon(self, i: int) -> Callable:
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
            >>> H: KHive = KHive(n=3, alpha=[3, 2, 0], beta=[2, 3, 0], gamma=[0, 0, 0], Uij=[[1, 0], [0]])
            >>> H_alpha.epsilon(i=1)(H=H)
            1
        """

        def _epsilon(H: KHive) -> int:
            return self._epsilon(i=i, k=i + 1)(H=H)

        return _epsilon

    def _epsilon(self, i: int, k: int) -> Callable:
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
            >>> H: KHive = KHive(n=3, alpha=[3, 2, 0], beta=[2, 3, 0], gamma=[0, 0, 0], Uij=[[1, 0], [0]])
            >>> H_alpha._epsilon(i=1, k=1)(H=H)
            0

            >>> H_alpha = KHives(n=3, alpha=[3, 2, 0])
            >>> H: KHive = KHive(n=3, alpha=[3, 2, 0], beta=[2, 3, 0], gamma=[0, 0, 0], Uij=[[1, 0], [0]])
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


@dataclass(frozen=True)
class FundamentalKHives(KHives):
    def __post_init__(self) -> None:
        self.is_fundamental_weight()

    def is_fundamental_weight(self) -> None:
        if not all(x in [0, 1] for x in self.alpha):
            raise ValueError("Invaid value! alpha is not a fundamental weight.")

    def phi(self, i: int) -> Callable:
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
            return max(self.weight(H=H)(i=i), 0)

        return _phi

    def f(self, i: int) -> Callable:
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
        """

        def _f(H: KHive) -> Union[KHive, None]:
            if self.phi(i=i)(H=H) == 0:
                return None
            else:
                i_as_index: int = i - 1

                beta: List[int] = H.beta.copy()
                gamma: List[int] = [0] * self.n
                Uij: List[List[int]] = deepcopy(H.Uij)

                beta[i_as_index] += -1
                beta[i_as_index + 1] += 1

                act_point_as_index: int
                try:
                    Uxi: List[int] = H.get_Uji()[i_as_index - 1]
                    act_point_as_index = [Uji > 0 for Uji in Uxi].index(True)
                except ValueError:
                    act_point_as_index = i_as_index

                i_as_index_in_Uij: int = i_as_index - act_point_as_index - 1
                ip1_as_index_in_Uij: int = i_as_index - act_point_as_index
                if 0 <= (i_as_index_in_Uij):
                    Uij[act_point_as_index][i_as_index_in_Uij] += -1  # U_{ki}
                Uij[act_point_as_index][ip1_as_index_in_Uij] += 1  # U_{k,i+1}

            return KHive(n=self.n, alpha=self.alpha, beta=beta, gamma=gamma, Uij=Uij)

        return _f

    def epsilon(self, i: int) -> Callable:
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
            return max(-self.weight(H=H)(i=i), 0)

        return _epsilon

    def e(self, i: int):
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
            >>> H = KHive(n=3, alpha=[1, 1, 0], beta=[1, 0, 1], gamma=[0, 0, 0], Uij=[[0, 0], [1]])
            >>> H_alpha.e(i=2)(H=H)
            KHive(n=3, alpha=[1, 1, 0], beta=[1, 1, 0], gamma=[0, 0, 0], Uij=[[0, 0], [0]])
        """

        def _e(H: KHive) -> Union[KHive, None]:
            if self.epsilon(i=i)(H=H) == 0:
                return None
            else:
                i_as_index: int = i - 1

                beta: List[int] = H.beta.copy()
                gamma: List[int] = [0] * self.n
                Uij: List[List[int]] = deepcopy(H.Uij)

                beta[i_as_index] += 1
                beta[i_as_index + 1] += -1

                act_point_as_index: int
                try:
                    Uxip1: List[int] = H.get_Uji()[i_as_index]
                    act_point_as_index = [Uji > 0 for Uji in Uxip1].index(True)
                except ValueError:
                    act_point_as_index = i_as_index

                i_as_index_in_Uij: int = i_as_index - act_point_as_index - 1
                ip1_as_index_in_Uij: int = i_as_index - act_point_as_index
                if 0 <= (i_as_index_in_Uij):
                    Uij[act_point_as_index][i_as_index_in_Uij] += +1  # U_{ki}
                Uij[act_point_as_index][ip1_as_index_in_Uij] += -1  # U_{k,i+1}

            return KHive(n=self.n, alpha=self.alpha, beta=beta, gamma=gamma, Uij=Uij)

        return _e
