import pytest

from pyloparsimony.examples import EXAMPLES
from pyloparsimony.parsimony import parsimony
from pyloparsimony.util import print_scenario, matrix_from_chars


def test_parsimony():
    ex1 = EXAMPLES["e1"]
    _ = parsimony(
        ex1['tree'],
        dict(zip(ex1['taxa'], ex1['patterns'])),
    )


def test_parsimony_up(capsys):
    ex1 = EXAMPLES["e1"]
    scenarios = parsimony(
        ex1['tree'],
        dict(zip(ex1['taxa'], ex1['patterns'])),
        characters=ex1['characters'],
        matrix=ex1['matrix'])

    print_scenario(scenarios[0], tree=ex1["tree"])
    out, _ = capsys.readouterr()
    assert 'Edge3/b' in out


@pytest.mark.parametrize(
    'chars,weights,kw,cond',
    [
        (
            ["a", "b", "c"],
            {("a", "b"): 2, ("a", "c"): 3},
            {},
            lambda r: r[0][1] != r[1][0]),
        (
            ["a", "b", "c"],
            {("a", "b"): 2, ("a", "c"): 3},
            dict(symmetric=True),
            lambda r: r[0][1] == r[1][0]),
    ]
)
def test_matrix_from_chars(chars, weights, kw, cond):
    assert cond(matrix_from_chars(chars, weights=weights, **kw))
