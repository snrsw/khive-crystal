import pytest
from khive_crystal.integer_hive_graph import IntegerHiveGraph


def test_fails_integer_hive_graph_invalid_uij():

    with pytest.raises(ValueError) as e:
        IntegerHiveGraph(
            n=3, alpha=[2, 1, 0], beta=[2, 1, 0], gamma=[0, 0, 0], uij=[[1, 0], [0]]
        )


def test_passes_integer_hive_graph_invalid_uij():

    assert isinstance(
        IntegerHiveGraph(
            n=3, alpha=[2, 1, 0], beta=[2, 1, 0], gamma=[0, 0, 0], uij=[[0, 0], [0]]
        ),
        IntegerHiveGraph,
    )
