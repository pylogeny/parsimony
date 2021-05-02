"""
Utility functions for parsimony calculations.
"""
from itertools import combinations
from pylotree import Tree

def matrix_from_chars(
        chars,
        model=None,
        weights=None,
        symmetric=False,
        default=1
        ):
    if model == "fitch":
        return [[1 if i != j else 0 for i in range(len(chars))] for j in range(len(chars))]

    if weights:
        matrix = [[0 for char in chars] for char in chars]
        if symmetric:
            for (i, charA), (j, charB) in combinations(enumerate(chars), r=2):
                matrix[i][j] = matrix[j][i] = weights.get(
                        (charA, charB),
                        weights.get((charB, charA), default))
        else:
            for i, charA in enumerate(chars):
                for j, charB in enumerate(chars):
                    if i != j:
                        matrix[i][j] = weights.get((charA, charB), default)
        return matrix


def print_scenario(scenario, tree):
    
    if isinstance(tree, Tree):
        tree = Tree(tree.newick)
    else:
        tree = Tree(tree)
        
    for edge, char in sorted(scenario, key=lambda x: len(x[0]), reverse=True):
        tree[edge].name = edge+'/'+char

    print(tree.ascii_art)
            

