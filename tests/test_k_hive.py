import pytest
from khive_crystal.k_hive import KHive


def test_fails_k_hive():
    with pytest.raises(ValueError) as e:
        KHive(n=3, alpha=[2, 1, 0], beta=[2, 1, 0], gamma=[0, 0, 0], Uij=[[1, 0], [0]])


def test_passes_k_hive():
    assert isinstance(
        KHive(n=3, alpha=[2, 1, 0], beta=[2, 1, 0], gamma=[0, 0, 0], Uij=[[0, 0], [0]]),
        KHive,
    )
