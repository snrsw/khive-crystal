from collections.abc import Generator
from typing import Any


def flatten_list(nestted_list: list[Any]) -> Generator[Any, None, None]:
    for element in nestted_list:
        if isinstance(element, list):
            yield from flatten_list(element)
        else:
            yield element
