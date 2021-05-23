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
        "matrix": matrix_from_chars(["a", "b", "c"])
    },
    "e2": {
        "characters": ["a", "b", "c"],
        "patterns": {
            "1": {
                "A": "b", 
                "B": "c", 
                "C": "a", 
                "D": "a", 
                "E": "b"
                },
            "2": {
                "A": ["b", "a"], 
                "B": "b", 
                "C": "c", 
                "D": "c", 
                "E": "b"
                }
            },
        "taxa": ["A", "B", "C", "D", "E"],
        "tree": "(((A,B)Edge1,(C,D)Edge2)Edge3,E)Root;",
        "matrix": matrix_from_chars(["a", "b", "c"])
        }
}
