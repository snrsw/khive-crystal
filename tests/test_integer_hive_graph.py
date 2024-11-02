import pytest

from khive_crystal.integer_hive_graph import IntegerHiveGraph


def test_fails_integer_hive_graph_invalid_Uij() -> None:

    with pytest.raises(ValueError):
        IntegerHiveGraph(
            n=3, alpha=[2, 1, 0], beta=[2, 1, 0], gamma=[0, 0, 0], Uij=[[1, 0], [0]]
        )


def test_passes_integer_hive_graph_invalid_Uij() -> None:

    assert isinstance(
        IntegerHiveGraph(
            n=3, alpha=[2, 1, 0], beta=[2, 1, 0], gamma=[0, 0, 0], Uij=[[0, 0], [0]]
        ),
        IntegerHiveGraph,
    )
