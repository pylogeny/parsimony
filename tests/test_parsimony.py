import pytest

from pyloparsimony.examples import EXAMPLES
from pyloparsimony.parsimony import parsimony, parsimony_analysis
from pyloparsimony.util import print_scenario, matrix_from_chars
from pylotree import Tree

def test_parsimony():
    ex1 = EXAMPLES["e1"]
    _ = parsimony(
        ex1['tree'],
        dict(zip(ex1['taxa'], ex1['patterns'])),
    )


def test_parsimony_analysis():
    ex2 = EXAMPLES["e2"]
    ex3 = EXAMPLES["e3"]
    assert parsimony_analysis(
            Tree(ex2["tree"]),
            ex2["patterns"]
            ) == 3
    assert parsimony_analysis(
            Tree(ex3["tree"]),
            ex3["patterns"],
            characters=ex3["characters"],
            matrices=ex3["matrices"]
            ) == 3



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


@pytest.mark.parametrize(
    'tree,pattern,root_states',
    [
        ('(A,B,C)', dict(A='a', B='a', C='b'), {'a'}),
        ('(A,B,C)', dict(A='a', B='b', C='c'), {'a', 'b', 'c'}),
        ('((A,B),C)', dict(A='a', B='a', C='b'), {'a', 'b'}),
    ]
)
def test_parsimony2(tree, pattern, root_states):
    assert set(dict(s)['Root'] for s in parsimony(tree, pattern)) == root_states
