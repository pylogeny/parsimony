"""
Utility functions for parsimony calculations.
"""
import itertools
import statistics
import random
from functools import partial

from pylotree import Tree

__all__ = ['matrix_from_chars', 'print_scenario', 'pointbiserialr']


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


def scenario_ascii_art(scenario, tree):
    tree = Tree(tree)
    scenario = dict(scenario)

    def relabel(n):
        if n.name in scenario:
            n.name = '{}/{}'.format(n.name, scenario[n.name])

    tree.root.visit(relabel)
    return tree.root.ascii_art()


def pvalues_by_iteration(x, y, function=None, iterate=1000):
    """
    Obtain significance scores for correlations functions through permutation tests.
    """
    r = function(x, y)
    vals = []
    xa = [a for a in x]
    for i in range(iterate):
        random.shuffle(xa)
        vals += [function(xa, y)]
    return r, len([v for v in vals if abs(v) > abs(r)]) / iterate


def pointbiserial_coefficient(x, y):
    """
    Compute the Point-Biserial Correlation Coefficient.

    .. note::
       
       Implementation follows the description in Wikipedia 
       (https://en.wikipedia.org/wiki/Point-biserial_correlation_coefficient).
    """

    assert len(x) == len(y)
    assert len(set(y)) == 2
    n = len(x)
    n1, n0 = y.count(1), y.count(0)
    m1 = statistics.mean([a for a, b in zip(x, y) if b == 1])
    m0 = statistics.mean([a for a, b in zip(x, y) if b == 0])
    xbar = statistics.mean(x)
    sn = (1/n * (sum([(xi - xbar) ** 2 for xi in x]))) ** 0.5
    r = ((m1-m0)/sn) * (((n1*n0)/(n**2))**0.5)
    return r


pointbiserialr = partial(pvalues_by_iteration,
        function=pointbiserial_coefficient)


