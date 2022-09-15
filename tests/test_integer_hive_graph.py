import pytest
from khive_crystal.integer_hive_graph import IntegerHiveGraph


def test_fails_integer_hive_graph_invalid_uij():

    with pytest.raises(ValueError) as e:
        IntegerHiveGraph(
            n=3, alpha=[1, 1, 1], beta=[1, 1, 1], gamma=[0, 0, 0], uij=[[1, 0], [0]]
        )

    assert (
        str(e.value)
        == "Invalid values! Given values dose not satisfy the definition of a interger hive graph."
    )


def test_passes_integer_hive_graph_invalid_uij():

    assert isinstance(
        IntegerHiveGraph(
            n=3, alpha=[1, 1, 1], beta=[1, 1, 1], gamma=[0, 0, 0], uij=[[0, 0], [0]]
        ),
        IntegerHiveGraph,
    )
