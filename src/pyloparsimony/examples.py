"""
Examples for parsimony calculation.
"""
from pyloparsimony.util import matrix_from_chars

EXAMPLES = {
        "e1": {
            "characters": ["a", "b", "c"],
            "patterns": ["b", "c", "a", "a", "b"],
            "taxa": ["A", "B", "C", "D", "E"],
            "tree": "(((A,B)Edge1,(C,D)Edge2)Edge3,E)Root;",
            "matrix": matrix_from_chars(["a", "b", "c"], model="fitch")
            },
        }
