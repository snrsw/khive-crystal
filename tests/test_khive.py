import pytest

from khive_crystal.khive import KHive


def test_fails_k_hive() -> None:
    with pytest.raises(ValueError):
        KHive(n=3, alpha=[1, 2, 0], beta=[1, 2, 0], gamma=[0, 0, 0], Uij=[[0, 0], [0]])

    with pytest.raises(ValueError):
        KHive(n=3, alpha=[3, 2, 0], beta=[3, 1, 0], gamma=[0, 1, 0], Uij=[[0, 0], [0]])

    with pytest.raises(ValueError):
        KHive(n=3, alpha=[3, 2, 0], beta=[3, 1, 0], gamma=[0, 0, 0], Uij=[[0, 0], [-1]])


def test_passes_k_hive() -> None:
    assert isinstance(
        KHive(n=3, alpha=[2, 1, 0], beta=[2, 1, 0], gamma=[0, 0, 0], Uij=[[0, 0], [0]]),
        KHive,
    )
