from typing import Generator


def flatten_list(nestted_list: list) -> Generator:
    for element in nestted_list:
        if isinstance(element, list):
            yield from flatten_list(element)
        else:
            yield element
