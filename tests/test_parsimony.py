from pyloparsimony.examples import EXAMPLES
from pyloparsimony.parsimony import parsimony, parsimony_up, parsimony_down

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


