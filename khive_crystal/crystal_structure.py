from collections.abc import Callable

from khive_crystal.fundamental_khives import FundamentalKHives
from khive_crystal.khive import KHive
from khive_crystal.khives import KHives
from khive_crystal.tensor_products import TensorProductsOfKHives


def phi(i: int) -> Callable[[KHive | list[KHive]], int]:
    """The entry point of phi_i. Given H and i, return phi_i(H) by appropriate crystal structure.

    Args:
        i (int): i in I

    Returns:
        Callable[[Union[KHive, List[KHive]]], int]: phi_i

    Examples:
        >>> H: KHive = KHive(n=3, alpha=[3, 2, 0], beta=[3, 2, 0], gamma=[0, 0, 0], Uij=[[0, 0], [0]])
        >>> phi(i=1)(H)
        1

        >>> H: KHive = KHive(n=3, alpha=[1, 1, 0], beta=[1, 1, 0], gamma=[0, 0, 0], Uij=[[0, 0], [0]])
        >>> phi(i=2)(H)
        1

        >>> H: KHive = KHive(n=3, alpha=[1, 1, 0], beta=[1, 1, 0], gamma=[0, 0, 0], Uij=[[0, 0], [0]])
        >>> phi(i=2)([H, H])
        2
    """  # noqa: B950

    def phi_i(H: KHive | list[KHive]) -> int:
        if isinstance(H, KHive):
            if H.is_fundamental_khive():
                return FundamentalKHives(n=H.n, alpha=H.alpha).phi(i=i)(H)
            else:
                return KHives(n=H.n, alpha=H.alpha).phi(i=i)(H)
        else:  # if H = List[KHives]
            tensor_products_khives: list[KHives | FundamentalKHives]
            if all(Hi.is_fundamental_khive() for Hi in H):
                tensor_products_khives = [
                    FundamentalKHives(n=Hi.n, alpha=Hi.alpha) for Hi in H
                ]
            else:
                tensor_products_khives = [KHives(n=Hi.n, alpha=Hi.alpha) for Hi in H]
            return TensorProductsOfKHives(
                tensor_products_khives=tensor_products_khives
            ).phi(i=i)(H)

    return phi_i


def f(
    i: int,
) -> Callable[
    [KHive | list[KHive]], KHive | None | list[KHive | None] | None
]:
    """The entry point of f_i. Given H and i, return f_i(H) by appropriate crystal structure.

    Args:
        i (int): i in I

    Returns:
        Callable[[Union[KHive, List[KHive]]], int]: f_i

    Examples:
        >>> H: KHive = KHive(n=3, alpha=[3, 2, 0], beta=[3, 2, 0], gamma=[0, 0, 0], Uij=[[0, 0], [0]])
        >>> f(i=1)(H)
        KHive(n=3, alpha=[3, 2, 0], beta=[2, 3, 0], gamma=[0, 0, 0], Uij=[[1, 0], [0]])

        >>> H: KHive = KHive(n=3, alpha=[1, 1, 0], beta=[1, 1, 0], gamma=[0, 0, 0], Uij=[[0, 0], [0]])
        >>> f(i=2)(H)
        KHive(n=3, alpha=[1, 1, 0], beta=[1, 0, 1], gamma=[0, 0, 0], Uij=[[0, 0], [1]])

        >>> H: KHive = KHive(n=3, alpha=[1, 1, 0], beta=[1, 1, 0], gamma=[0, 0, 0], Uij=[[0, 0], [0]])
        >>> f(i=2)([H, H])
        [KHive(n=3, alpha=[1, 1, 0], beta=[1, 0, 1], gamma=[0, 0, 0], Uij=[[0, 0], [1]]), KHive(n=3, alpha=[1, 1, 0], beta=[1, 1, 0], gamma=[0, 0, 0], Uij=[[0, 0], [0]])]
    """  # noqa: E501, B950

    def f_i(
        H: KHive | list[KHive]
    ) -> KHive | None | list[KHive | None] | None:
        if isinstance(H, KHive):
            if H.is_fundamental_khive():
                return FundamentalKHives(n=H.n, alpha=H.alpha).f(i=i)(H)
            else:
                return KHives(n=H.n, alpha=H.alpha).f(i=i)(H)
        else:  # if H = List[KHives]
            tensor_products_khives: list[KHives | FundamentalKHives]
            if all(Hi.is_fundamental_khive() for Hi in H):
                tensor_products_khives = [
                    FundamentalKHives(n=Hi.n, alpha=Hi.alpha) for Hi in H
                ]
            else:
                tensor_products_khives = [KHives(n=Hi.n, alpha=Hi.alpha) for Hi in H]
            return TensorProductsOfKHives(
                tensor_products_khives=tensor_products_khives
            ).f(i=i)(H)

    return f_i


def epsilon(i: int) -> Callable[[KHive | list[KHive]], int]:
    """The entry point of epsilon_i.
    Given H and i, return epsilon_i(H) by appropriate crystal structure.

    Args:
        i (int): i in I

    Returns:
        Callable[[Union[KHive, List[KHive]]], int]: epsilon_i

    Examples:
        >>> H: KHive = KHive(
        ...    n=3,
        ...    alpha=[3, 2, 0],
        ...    beta=[2, 3, 0],
        ...    gamma=[0, 0, 0],
        ...    Uij=[[1, 0], [0]]
        ... )
        >>> epsilon(i=1)(H)
        1

        >>> H: KHive = KHive(n=3, alpha=[1, 1, 0], beta=[1, 1, 0], gamma=[0, 0, 0], Uij=[[0, 0], [0]])
        >>> epsilon(i=2)(H)
        0

        >>> H: List[KHive] = [
        ...    KHive(n=3, alpha=[1, 1, 0], beta=[1, 0, 1], gamma=[0, 0, 0], Uij=[[0, 0], [1]]),
        ...    KHive(n=3, alpha=[1, 1, 0], beta=[1, 1, 0], gamma=[0, 0, 0], Uij=[[0, 0], [0]]),
        ... ]
        >>> epsilon(i=2)(H)
        1
    """  # noqa: B950

    def epsilon_i(H: KHive | list[KHive]) -> int:
        if isinstance(H, KHive):
            if H.is_fundamental_khive():
                return FundamentalKHives(n=H.n, alpha=H.alpha).epsilon(i=i)(H)
            else:
                return KHives(n=H.n, alpha=H.alpha).epsilon(i=i)(H)
        else:  # if H = List[KHives]
            tensor_products_khives: list[KHives | FundamentalKHives]
            if all(Hi.is_fundamental_khive() for Hi in H):
                tensor_products_khives = [
                    FundamentalKHives(n=Hi.n, alpha=Hi.alpha) for Hi in H
                ]
            else:
                tensor_products_khives = [KHives(n=Hi.n, alpha=Hi.alpha) for Hi in H]
            return TensorProductsOfKHives(
                tensor_products_khives=tensor_products_khives
            ).epsilon(i=i)(H)

    return epsilon_i


def e(
    i: int,
) -> Callable[
    [KHive | list[KHive]], KHive | None | list[KHive | None] | None
]:
    """The entry point of e_i. Given H and i, return e_i(H) by appropriate crystal structure.

    Args:
        i (int): i in I

    Returns:
        Callable[[Union[KHive, List[KHive]]], int]: f_i

    Examples:
        >>> H: KHive = KHive(
        ...    n=3,
        ...    alpha=[3, 2, 0],
        ...    beta=[2, 3, 0],
        ...    gamma=[0, 0, 0],
        ...    Uij=[[1, 0], [0]]
        ... )
        >>> e(i=1)(H)
        KHive(n=3, alpha=[3, 2, 0], beta=[3, 2, 0], gamma=[0, 0, 0], Uij=[[0, 0], [0]])

        >>> H = KHive(
        ...     n=3,
        ...     alpha=[1, 1, 0],
        ...     beta=[1, 0, 1],
        ...     gamma=[0, 0, 0],
        ...     Uij=[[0, 0], [1]]
        ... )
        >>> e(i=2)(H)
        KHive(n=3, alpha=[1, 1, 0], beta=[1, 1, 0], gamma=[0, 0, 0], Uij=[[0, 0], [0]])

        >>> H: KHive = KHive(n=3, alpha=[1, 1, 0], beta=[1, 1, 0], gamma=[0, 0, 0], Uij=[[0, 0], [0]])
        >>> e(i=2)([H, H])

    """  # noqa: B950

    def e_i(
        H: KHive | list[KHive]
    ) -> KHive | None | list[KHive | None] | None:
        if isinstance(H, KHive):
            if H.is_fundamental_khive():
                return FundamentalKHives(n=H.n, alpha=H.alpha).e(i=i)(H)
            else:
                return KHives(n=H.n, alpha=H.alpha).e(i=i)(H)
        else:  # if H = List[KHives]
            tensor_products_khives: list[KHives | FundamentalKHives]
            if all(Hi.is_fundamental_khive() for Hi in H):
                tensor_products_khives = [
                    FundamentalKHives(n=Hi.n, alpha=Hi.alpha) for Hi in H
                ]
            else:
                tensor_products_khives = [KHives(n=Hi.n, alpha=Hi.alpha) for Hi in H]
            return TensorProductsOfKHives(
                tensor_products_khives=tensor_products_khives
            ).e(i=i)(H)

    return e_i
