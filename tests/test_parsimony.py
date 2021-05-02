from pyloparsimony.examples import EXAMPLES
from pyloparsimony.parsimony import parsimony, parsimony_up, parsimony_down
from pyloparsimony.util import print_scenario, matrix_from_chars
from pylotree.tree import Tree

def test_parsimony():

    ex1 = EXAMPLES["e1"]
    weights, scenarios, weight_dict = parsimony(**ex1)

def test_parsimony_up():
    ex1 = EXAMPLES["e1"]
    weight_dict = parsimony_up(**ex1)
    ex1["output"] = "weight"
    weight = parsimony_up(**ex1)
    assert weight == 2
    ex1["output"] = "chars"
    chars = parsimony_up(**ex1)

    scenarios = parsimony_down(
            weights=weight_dict,
            patterns=ex1["patterns"],
            taxa=ex1["taxa"],
            tree=ex1["tree"],
            matrix=ex1["matrix"],
            characters=ex1["characters"],
            )
    print_scenario(scenarios[0], tree=ex1["tree"])
    print_scenario(scenarios[0], tree=Tree(ex1["tree"]))
    

def test_matrix_from_chars():

    matrix = matrix_from_chars(["a", "b", "c"], weights={
        ("a", "b"): 2,
        ("a", "c"): 3
        })
    assert matrix[0][1] == 2
    assert matrix[0][2] == 3
    assert matrix[1][0] == 1

    matrix = matrix_from_chars(["a", "b", "c"], weights={
        ("a", "b"): 2,
        ("a", "c"): 3
        },
        symmetric=True)
    assert matrix[0][1] == 2
    assert matrix[0][2] == 3
    assert matrix[1][0] == 2

    






