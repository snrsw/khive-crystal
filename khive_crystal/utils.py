from typing import Generator, List


def flatten_list(nestted_list: list) -> Generator:
    for element in nestted_list:
        if isinstance(element, list):
            yield from flatten_list(element)
        else:
            yield element


def get_length(alpha: List[int]) -> int:
    """_summary_

    Args:
        alpha (List[int]): A composition.

    Returns:
        int:  length of alpha.

    Examples:
    >>> get_length(alpha=[3, 2, 1, 0])
    3
    >>> get_length(alpha=[3, 2, 1, 0, 0])
    3
    >>> get_length(alpha=[3, 2, 0, 0, 0])
    2
    """
    return [alpha_i > 0 for alpha_i in alpha].index(False)
