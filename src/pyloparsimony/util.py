"""
Utility functions for parsimony calculations.
"""
import itertools

from pylotree import Tree

__all__ = ['matrix_from_chars', 'print_scenario']


def matrix_from_chars(chars, weights=None, symmetric=False, default=1):
    """
    With the default argument values, this returns a matrix according to the "fitch" model.
    """
    weights = weights or {}
    matrix = [[0 for _ in chars] for _ in chars]
    if symmetric:
        for (i, charA), (j, charB) in itertools.combinations(enumerate(chars), r=2):
            matrix[i][j] = matrix[j][i] = weights.get(
                    (charA, charB),
                    weights.get((charB, charA), default))
    else:
        for (i, charA), (j, charB) in itertools.permutations(enumerate(chars), 2):
            matrix[i][j] = weights.get((charA, charB), default)
    return matrix


def print_scenario(scenario, tree):
    tree = Tree(tree)
    scenario = dict(scenario)

    def relabel(n):
        if n.name in scenario:
            n.name = '{}/{}'.format(n.name, scenario[n.name])

    tree.root.visit(relabel)
    print(tree.root.ascii_art())




