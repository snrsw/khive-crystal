from logging import DEBUG, StreamHandler, getLogger
from typing import Callable, List, Union

from khive_crystal.khive import KHive
from khive_crystal.khives import KHives
from khive_crystal.utils import flatten_list

logger = getLogger(__name__)
handler = StreamHandler()
handler.setLevel(DEBUG)
logger.setLevel(DEBUG)
logger.addHandler(handler)
logger.propagate = False


class TensorProductsOfKHives:
    def __init__(self, tensort_products_khives: List[KHives]) -> None:
        self.tensort_products_khives: List[KHives] = tensort_products_khives

    def phi(self, i: int) -> Callable[[List[KHive]], int]:
        """phi_i(tensort_products_khive)

        Args:
            i (int): i in I

        Returns:
            Callable[[List[KHive]], int]: phi_i

        Examples:
            >>> khives: KHives = KHives(n=3, alpha=[1, 1, 0])
            >>> tensor: TensorProductsOfKHives = TensorProductsOfKHives(
            ...    tensort_products_khives=[khives, khives]
            ... )
            >>> H: List[KHive] = [
            ...    khives.highest_weight_vector(),
            ...    khives.highest_weight_vector()
            ... ]
            >>> tensor.phi(i=2)(H)
            2

            >>> khives: KHives = KHives(n=3, alpha=[1, 1, 0])
            >>> tensor: TensorProductsOfKHives = TensorProductsOfKHives(
            ...    tensort_products_khives=[khives, khives]
            ... )
            >>> H: List[KHive] = [
            ...    khives.highest_weight_vector(),
            ...    khives.highest_weight_vector()
            ... ]
            >>> tensor.phi(i=1)(H)
            0
        """

        def phi_i(tensort_products_khive: List[KHive]) -> int:
            left_khives: KHives = self.tensort_products_khives[0]
            left_khive: KHive = tensort_products_khive[0]

            right_khives: Union[TensorProductsOfKHives, KHives]
            right_khive: Union[List[KHive], KHive]
            if len(self.tensort_products_khives[1:]) == 1:
                right_khives = self.tensort_products_khives[1]
                right_khive = tensort_products_khive[1]

                return max(
                    right_khives.phi(i=i)(right_khive),
                    left_khives.phi(i=i)(left_khive)
                    + right_khives.inner_product(
                        i=i, weight=right_khives.weight(H=right_khive)
                    ),
                )

            right_khives = TensorProductsOfKHives(
                tensort_products_khives=self.tensort_products_khives[1:]
            )
            right_khive = tensort_products_khive[1:]

            return max(
                right_khives.phi(i=i)(right_khive),
                left_khives.phi(i=i)(left_khive)
                + right_khives.inner_product(
                    i=i, weight=right_khives.weight(tensort_products_khive=right_khive)
                ),
            )

        return phi_i

    def f(self, i: int) -> Callable[[List[KHive]], List[KHive]]:
        """f_i(tensort_products_khive)

        Args:
            i (int): i in I

        Returns:
            Callable[[List[KHive]], List[KHive]]: f_i

        Examples:
            >>> khives: KHives = KHives(n=3, alpha=[1, 1, 0])
            >>> tensor: TensorProductsOfKHives = TensorProductsOfKHives(
            ...    tensort_products_khives=[khives, khives]
            ... )
            >>> H: List[KHive] = [
            ...    khives.highest_weight_vector(),
            ...    khives.highest_weight_vector()
            ... ]
            >>> tensor.f(i=2)(H)
            [KHive(n=3, alpha=[1, 1, 0], beta=[1, 0, 1], gamma=[0, 0, 0], Uij=[[0, 0], [1]]), KHive(n=3, alpha=[1, 1, 0], beta=[1, 1, 0], gamma=[0, 0, 0], Uij=[[0, 0], [0]])]
        """  # noqa: B950

        def f_i(tensort_products_khive: List[KHive]) -> List[KHive]:
            left_khives: KHives = self.tensort_products_khives[0]
            left_khive: KHive = tensort_products_khive[0]

            right_khives: Union[TensorProductsOfKHives, KHives]
            right_khive: Union[List[KHive], KHive]
            right_epsilon_i: int
            left_phi_i: int
            if len(self.tensort_products_khives[1:]) == 1:
                right_khives = self.tensort_products_khives[1]
                right_khive = tensort_products_khive[1]
                left_phi_i = left_khives.phi(i=i)(left_khive)
                right_epsilon_i = right_khives.epsilon(i=i)(right_khive)
                if left_phi_i > right_epsilon_i:
                    return list(
                        flatten_list([left_khives.f(i=i)(left_khive), right_khive])
                    )
                else:
                    return list(
                        flatten_list([left_khive, right_khives.f(i=1)(right_khive)])
                    )

            right_khives = TensorProductsOfKHives(
                tensort_products_khives=self.tensort_products_khives[1:]
            )
            right_khive = tensort_products_khive[1:]
            left_phi_i = left_khives.phi(i=i)(left_khive)
            right_epsilon_i = right_khives.epsilon(i=i)(right_khive)

            if left_phi_i > right_epsilon_i:
                return list(flatten_list([left_khives.f(i=i)(left_khive), right_khive]))
            else:
                return list(
                    flatten_list([left_khive, right_khives.f(i=1)(right_khive)])
                )

        return f_i

    def epsilon(self, i: int) -> Callable[[List[KHive]], int]:
        """epsilon_i(tensort_products_khive)

        Args:
            i (int): i in I

        Returns:
            Callable[[List[KHive]], int]: epsilon_i

        Examples:
            >>> khives: KHives = KHives(n=3, alpha=[1, 1, 0])
            >>> tensor: TensorProductsOfKHives = TensorProductsOfKHives(
            ...    tensort_products_khives=[khives, khives]
            ... )
            >>> H: List[KHive] = [
            ...    khives.highest_weight_vector(),
            ...    khives.highest_weight_vector()
            ... ]
            >>> tensor.epsilon(i=2)(H)
            0

            >>> khives: KHives = KHives(n=3, alpha=[1, 1, 0])
            >>> tensor: TensorProductsOfKHives = TensorProductsOfKHives(
            ...    tensort_products_khives=[khives, khives]
            ... )
            >>> H: List[KHive] = [
            ...    KHive(n=3, alpha=[1, 1, 0], beta=[1, 0, 1], gamma=[0, 0, 0], Uij=[[0, 0], [1]]),
            ...    khives.highest_weight_vector()
            ... ]
            >>> tensor.epsilon(i=2)(H)
            1
        """  # noqa: B950

        def epsilon_i(tensort_products_khive: List[KHive]) -> int:
            left_khives: KHives = self.tensort_products_khives[0]
            left_khive: KHive = tensort_products_khive[0]

            right_khives: Union[TensorProductsOfKHives, KHives]
            right_khive: Union[List[KHive], KHive]
            if len(self.tensort_products_khives[1:]) == 1:
                right_khives = self.tensort_products_khives[1]
                right_khive = tensort_products_khive[1]

                return max(
                    left_khives.epsilon(i=i)(left_khive),
                    right_khives.epsilon(i=i)(right_khive)
                    - left_khives.inner_product(
                        i=i, weight=left_khives.weight(H=left_khive)
                    ),
                )

            right_khives = TensorProductsOfKHives(
                tensort_products_khives=self.tensort_products_khives[1:]
            )
            right_khive = tensort_products_khive[1:]

            return max(
                left_khives.epsilon(i=i)(left_khive),
                right_khives.epsilon(i=i)(right_khive)
                - left_khives.inner_product(
                    i=i,
                    weight=left_khives.weight(H=left_khive),
                ),
            )

        return epsilon_i

    def e(self, i: int) -> Callable[[List[KHive]], List[KHive]]:
        """f_i(tensort_products_khive)

        Args:
            i (int): i in I

        Returns:
            Callable[[List[KHive]], List[KHive]]: f_i

        Examples:
            >>> khives: KHives = KHives(n=3, alpha=[1, 1, 0])
            >>> tensor: TensorProductsOfKHives = TensorProductsOfKHives(
            ...    tensort_products_khives=[khives, khives]
            ... )
            >>> H: List[KHive] = [
            ...    khives.highest_weight_vector(),
            ...    khives.highest_weight_vector()
            ... ]
            >>> tensor.e(i=2)(H)
            [None, None]
        """  # noqa: B950

        def e_i(tensort_products_khive: List[KHive]) -> List[KHive]:
            left_khives: KHives = self.tensort_products_khives[0]
            left_khive: KHive = tensort_products_khive[0]

            right_khives: Union[TensorProductsOfKHives, KHives]
            right_khive: Union[List[KHive], KHive]
            right_epsilon_i: int
            left_phi_i: int
            if len(self.tensort_products_khives[1:]) == 1:
                right_khives = self.tensort_products_khives[1]
                right_khive = tensort_products_khive[1]
                left_phi_i = left_khives.phi(i=i)(left_khive)
                right_epsilon_i = right_khives.epsilon(i=i)(right_khive)
                if left_phi_i >= right_epsilon_i:
                    return list(
                        flatten_list([left_khives.e(i=i)(left_khive), right_khive])
                    )
                else:
                    return list(
                        flatten_list([left_khive, right_khives.e(i=1)(right_khive)])
                    )

            right_khives = TensorProductsOfKHives(
                tensort_products_khives=self.tensort_products_khives[1:]
            )
            right_khive = tensort_products_khive[1:]
            left_phi_i = left_khives.phi(i=i)(left_khive)
            right_epsilon_i = right_khives.epsilon(i=i)(right_khive)

            if left_phi_i >= right_epsilon_i:
                return list(flatten_list([left_khives.e(i=i)(left_khive), right_khive]))
            else:
                return list(
                    flatten_list([left_khive, right_khives.e(i=1)(right_khive)])
                )

        return e_i

    def weight(self, tensort_products_khive: List[KHive]) -> List[int]:
        """weight(tensort_products_khive)

        Args:
            tensort_products_khive (List[KHive]): _description_

        Returns:
            List[int]: weight(tensort_products_khive)

        Examples:
            >>> khives: KHives = KHives(n=3, alpha=[1, 1, 0])
            >>> tensor: TensorProductsOfKHives = TensorProductsOfKHives(
            ...    tensort_products_khives=[khives, khives]
            ... )
            >>> H: List[KHive] = [
            ...    khives.highest_weight_vector(),
            ...    khives.highest_weight_vector()
            ... ]
            >>> tensor.weight(H)
            [2, 2, 0]
        """
        weights: List[List[int]] = [
            khives.weight(H=H)
            for H, khives in zip(tensort_products_khive, self.tensort_products_khives)
        ]
        return [sum(weights_i) for weights_i in zip(*weights)]

    def inner_product(self, i: int, weight: List[int]) -> int:
        pass
