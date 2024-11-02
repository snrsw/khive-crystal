from collections.abc import Callable

from khive_crystal.khive import KHive
from khive_crystal.khives import KHives
from khive_crystal.utils import flatten_list


class TensorProductsOfKHives:
    """
    This class implements the tensor product rule for K-hives.
    """

    def __init__(self, tensor_products_khives: list[KHives]) -> None:
        self.tensor_products_khives: list[KHives] = tensor_products_khives

    def phi(self, i: int) -> Callable[[list[KHive]], int]:
        """phi_i(tensor_products_khive)

        Args:
            i (int): i in I

        Returns:
            Callable[[List[KHive]], int]: phi_i

        Examples:
            >>> khives: KHives = KHives(n=3, alpha=[1, 1, 0])
            >>> tensor: TensorProductsOfKHives = TensorProductsOfKHives(
            ...    tensor_products_khives=[khives, khives]
            ... )
            >>> H: List[KHive] = [
            ...    khives.highest_weight_vector(),
            ...    khives.highest_weight_vector()
            ... ]
            >>> tensor.phi(i=2)(H)
            2

            >>> khives: KHives = KHives(n=3, alpha=[1, 1, 0])
            >>> tensor: TensorProductsOfKHives = TensorProductsOfKHives(
            ...    tensor_products_khives=[khives, khives]
            ... )
            >>> H: List[KHive] = [
            ...    khives.highest_weight_vector(),
            ...    khives.highest_weight_vector()
            ... ]
            >>> tensor.phi(i=1)(H)
            0
        """

        def phi_i(tensor_products_khive: list[KHive]) -> int:
            left_khives: KHives = self.tensor_products_khives[0]
            left_khive: KHive = tensor_products_khive[0]

            right_khives: TensorProductsOfKHives | KHives
            right_khive: list[KHive] | KHive
            if len(self.tensor_products_khives[1:]) == 1:
                right_khives = self.tensor_products_khives[1]
                right_khive = tensor_products_khive[1]

                return max(
                    right_khives.phi(i=i)(right_khive),
                    left_khives.phi(i=i)(left_khive)
                    + right_khives.inner_product(
                        i=i, weight=right_khives.weight(H=right_khive)
                    ),
                )

            right_khives = TensorProductsOfKHives(
                tensor_products_khives=self.tensor_products_khives[1:]
            )
            right_khive = tensor_products_khive[1:]

            return max(
                right_khives.phi(i=i)(right_khive),
                left_khives.phi(i=i)(left_khive)
                + right_khives.inner_product(
                    i=i, weight=right_khives.weight(tensor_products_khive=right_khive)
                ),
            )

        return phi_i

    def f(self, i: int) -> Callable[[list[KHive]], list[KHive | None] | None]:
        """f_i(tensor_products_khive)

        Args:
            i (int): i in I

        Returns:
            Callable[[List[KHive]], List[KHive]]: f_i

        Examples:
            >>> khives: KHives = KHives(n=3, alpha=[1, 1, 0])
            >>> tensor: TensorProductsOfKHives = TensorProductsOfKHives(
            ...    tensor_products_khives=[khives, khives]
            ... )
            >>> H: List[KHive] = [
            ...    khives.highest_weight_vector(),
            ...    khives.highest_weight_vector()
            ... ]
            >>> tensor.f(i=2)(H)
            [KHive(n=3, alpha=[1, 1, 0], beta=[1, 0, 1], gamma=[0, 0, 0], Uij=[[0, 0], [1]]), KHive(n=3, alpha=[1, 1, 0], beta=[1, 1, 0], gamma=[0, 0, 0], Uij=[[0, 0], [0]])]
        """  # noqa: B950

        def f_i(tensor_products_khive: list[KHive]) -> list[KHive | None] | None:
            left_khives: KHives = self.tensor_products_khives[0]
            left_khive: KHive = tensor_products_khive[0]

            right_khives: TensorProductsOfKHives | KHives
            right_khive: list[KHive] | KHive
            right_epsilon_i: int
            left_phi_i: int
            if len(self.tensor_products_khives[1:]) == 1:
                right_khives = self.tensor_products_khives[1]
                right_khive = tensor_products_khive[1]
                left_phi_i = left_khives.phi(i=i)(left_khive)
                right_epsilon_i = right_khives.epsilon(i=i)(right_khive)
                if left_phi_i > right_epsilon_i:
                    return self.zero_scalar(
                        list(
                            flatten_list([left_khives.f(i=i)(left_khive), right_khive])
                        )
                    )
                else:
                    return self.zero_scalar(
                        list(
                            flatten_list([left_khive, right_khives.f(i=1)(right_khive)])
                        )
                    )

            right_khives = TensorProductsOfKHives(
                tensor_products_khives=self.tensor_products_khives[1:]
            )
            right_khive = tensor_products_khive[1:]
            left_phi_i = left_khives.phi(i=i)(left_khive)
            right_epsilon_i = right_khives.epsilon(i=i)(right_khive)

            if left_phi_i > right_epsilon_i:
                return self.zero_scalar(
                    list(flatten_list([left_khives.f(i=i)(left_khive), right_khive]))
                )
            else:
                return self.zero_scalar(
                    list(flatten_list([left_khive, right_khives.f(i=1)(right_khive)]))
                )

        return f_i

    def epsilon(self, i: int) -> Callable[[list[KHive]], int]:
        """epsilon_i(tensor_products_khive)

        Args:
            i (int): i in I

        Returns:
            Callable[[List[KHive]], int]: epsilon_i

        Examples:
            >>> khives: KHives = KHives(n=3, alpha=[1, 1, 0])
            >>> tensor: TensorProductsOfKHives = TensorProductsOfKHives(
            ...    tensor_products_khives=[khives, khives]
            ... )
            >>> H: List[KHive] = [
            ...    khives.highest_weight_vector(),
            ...    khives.highest_weight_vector()
            ... ]
            >>> tensor.epsilon(i=2)(H)
            0

            >>> khives: KHives = KHives(n=3, alpha=[1, 1, 0])
            >>> tensor: TensorProductsOfKHives = TensorProductsOfKHives(
            ...    tensor_products_khives=[khives, khives]
            ... )
            >>> H: List[KHive] = [
            ...    KHive(n=3, alpha=[1, 1, 0], beta=[1, 0, 1], gamma=[0, 0, 0], Uij=[[0, 0], [1]]),
            ...    khives.highest_weight_vector()
            ... ]
            >>> tensor.epsilon(i=2)(H)
            1
        """  # noqa: B950

        def epsilon_i(tensor_products_khive: list[KHive]) -> int:
            left_khives: KHives = self.tensor_products_khives[0]
            left_khive: KHive = tensor_products_khive[0]

            right_khives: TensorProductsOfKHives | KHives
            right_khive: list[KHive] | KHive
            if len(self.tensor_products_khives[1:]) == 1:
                right_khives = self.tensor_products_khives[1]
                right_khive = tensor_products_khive[1]

                return max(
                    left_khives.epsilon(i=i)(left_khive),
                    right_khives.epsilon(i=i)(right_khive)
                    - left_khives.inner_product(
                        i=i, weight=left_khives.weight(H=left_khive)
                    ),
                )

            right_khives = TensorProductsOfKHives(
                tensor_products_khives=self.tensor_products_khives[1:]
            )
            right_khive = tensor_products_khive[1:]

            return max(
                left_khives.epsilon(i=i)(left_khive),
                right_khives.epsilon(i=i)(right_khive)
                - left_khives.inner_product(
                    i=i,
                    weight=left_khives.weight(H=left_khive),
                ),
            )

        return epsilon_i

    def e(self, i: int) -> Callable[[list[KHive]], list[KHive | None] | None]:
        """f_i(tensor_products_khive)

        Args:
            i (int): i in I

        Returns:
            Callable[[List[KHive]], List[KHive]]: f_i

        Examples:
            >>> khives: KHives = KHives(n=3, alpha=[1, 1, 0])
            >>> tensor: TensorProductsOfKHives = TensorProductsOfKHives(
            ...    tensor_products_khives=[khives, khives]
            ... )
            >>> H: List[KHive] = [
            ...    khives.highest_weight_vector(),
            ...    khives.highest_weight_vector()
            ... ]
            >>> tensor.e(i=2)(H)
        """  # noqa: B950

        def e_i(tensor_products_khive: list[KHive]) -> list[KHive | None] | None:
            left_khives: KHives = self.tensor_products_khives[0]
            left_khive: KHive = tensor_products_khive[0]

            right_khives: TensorProductsOfKHives | KHives
            right_khive: list[KHive] | KHive
            right_epsilon_i: int
            left_phi_i: int
            if len(self.tensor_products_khives[1:]) == 1:
                right_khives = self.tensor_products_khives[1]
                right_khive = tensor_products_khive[1]
                left_phi_i = left_khives.phi(i=i)(left_khive)
                right_epsilon_i = right_khives.epsilon(i=i)(right_khive)
                if left_phi_i >= right_epsilon_i:
                    return self.zero_scalar(
                        list(
                            flatten_list([left_khives.e(i=i)(left_khive), right_khive])
                        )
                    )
                else:
                    return self.zero_scalar(
                        list(
                            flatten_list([left_khive, right_khives.e(i=1)(right_khive)])
                        )
                    )

            right_khives = TensorProductsOfKHives(
                tensor_products_khives=self.tensor_products_khives[1:]
            )
            right_khive = tensor_products_khive[1:]
            left_phi_i = left_khives.phi(i=i)(left_khive)
            right_epsilon_i = right_khives.epsilon(i=i)(right_khive)

            if left_phi_i >= right_epsilon_i:
                return self.zero_scalar(
                    list(flatten_list([left_khives.e(i=i)(left_khive), right_khive]))
                )
            else:
                return self.zero_scalar(
                    list(flatten_list([left_khive, right_khives.e(i=1)(right_khive)]))
                )

        return e_i

    def weight(self, tensor_products_khive: list[KHive]) -> list[int]:
        """weight(tensor_products_khive)

        Args:
            tensor_products_khive (List[KHive]): _description_

        Returns:
            List[int]: weight(tensor_products_khive)

        Examples:
            >>> khives: KHives = KHives(n=3, alpha=[1, 1, 0])
            >>> tensor: TensorProductsOfKHives = TensorProductsOfKHives(
            ...    tensor_products_khives=[khives, khives]
            ... )
            >>> H: List[KHive] = [
            ...    khives.highest_weight_vector(),
            ...    khives.highest_weight_vector()
            ... ]
            >>> tensor.weight(H)
            [2, 2, 0]
        """
        weights: list[list[int]] = [
            khives.weight(H=H)
            for H, khives in zip(tensor_products_khive, self.tensor_products_khives, strict=False)
        ]
        return [sum(weights_i) for weights_i in zip(*weights, strict=False)]

    def inner_product(self, i: int, weight: list[int]) -> int:
        """Compute <h_i, wt(H)>

        Args:
            i (int): h_i
            weight (List[int]): weight.

        Returns:
            int: <h_i, wt(H)>
        """
        if i not in [_ + 1 for _ in range(self.tensor_products_khives[0].n - 1)]:
            raise ValueError("i must be in I.")

        i_as_index: int = i - 1
        return weight[i_as_index] - weight[i_as_index + 1]

    def zero_scalar(
        self, tensor_products_khive: list[KHive | None]
    ) -> list[KHive | None] | None:
        """zero scalar

        Args:
            tensor_products_khive (List[Optional[KHive]]): Result of the tensor product rule.

        Returns:
            Optional[List[Optional[KHive]]]: If tensor_products_khive has None, return None

        Examples:
            >>> khives: KHives = KHives(n=3, alpha=[1, 1, 0])
            >>> tensor: TensorProductsOfKHives = TensorProductsOfKHives(
            ...    tensor_products_khives=[khives, khives]
            ... )
            >>> H: List[KHive] = [
            ...    None,
            ...    khives.highest_weight_vector()
            ... ]
            >>> tensor.zero_scalar(H)
        """
        return None if None in tensor_products_khive else tensor_products_khive
