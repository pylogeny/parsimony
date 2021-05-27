"""
Compute parsimony analyses
"""
import itertools
import collections

from pylotree import Tree
from .util import matrix_from_chars


def parsimony(tree, pattern, characters=None, matrix=None):
    """
    Calculate the most parsimonious evolutionary scenario for pattern and tree. 
    """
    tree = Tree(tree)
    # Take the set of all observed states as default for all possible characters:
    if not characters:
        characters = []
        for pt in pattern.values():
            if isinstance(pt, (list, tuple)):
                characters += pt
            else:
                characters += [pt]
        characters = sorted(set(characters))

    # Use a Fitch-model matrix as default:
    matrix = matrix or matrix_from_chars(characters)

    assert set(pattern.keys()).issubset(tree.root.get_leaf_names())

    return down(tree, characters, matrix, up(tree, characters, matrix, pattern))


def up(tree, characters, matrix, pattern):
    W = collections.defaultdict(dict)

    for node, name in map(lambda x: (x, x.name), tree.postorder):        
        if not node.descendants:
            for i, char in enumerate(characters):
                W[name][char] = 0 if char in pattern[name] else 1000000
        else:
            for i, charA in enumerate(characters):
                scoresA = []
                for nodeB, nameB in map(lambda x: (x, x.name), node.descendants):
                    scoresB = []
                    for j, charB in enumerate(characters):
                        weightB = W[nameB][charB]
                        weightNew = weightB + matrix[i][j]
                        scoresB += [weightNew]
                    smin = min(scoresB)
                    scoresA += [smin]
                W[name][charA] = sum(scoresA)

    return W


def down(tree, characters, matrix, weights):
    smin = min(weights[tree.name].values())
    root_chars = [a for a, b in weights[tree.name].items() if b == smin]

    descendants = {n.name: [nn for nn in n.descendants] for n in tree}

    # prepare the queue
    queue = []
    for char in root_chars:
        nodes = descendants[tree.name]
        queue += [([(nodes, tree.name, char, characters.index(char))], [(tree.name, char)])]
    
    output = []
    while queue:
        nodes, scenario = queue.pop(0)
        if not nodes:
            output += [scenario]
        else:
            children, parent, pchar, pidx = nodes.pop()
            pscore = weights[parent][pchar]
            for comb in itertools.product(*len(children) * [characters]):
                score = 0
                for i, char in enumerate(comb):
                    cidx = characters.index(char)
                    score += matrix[pidx][cidx]
                    score += weights[children[i].name][char]
                if score == pscore:
                    new_nodes = [n for n in nodes]
                    new_scenario = [s for s in scenario]
                    for child, char in zip(children,comb):
                        new_nodes += [
                            (descendants[child.name], child.name, char, characters.index(char))]
                        new_scenario += [(child.name, char)]
                    queue += [(new_nodes, new_scenario)]
    return output


def parsimony_analysis(tree, patterns, characters=None, matrices=None):
    """
    Carry out a parsimony analysis for a given tree and a number of patterns.
    """
    characters = characters or {}
    matrices = matrices or {}
    scenarios, weights, root_weights = {}, {}, {}
    for key, pattern in patterns.items():
        characterlist, matrix = characters.get(key), matrices.get(key)
        if not characterlist:
            chars = []
            for sublist in pattern.values():
                chars += sublist
            characterlist = sorted(set(chars))
        if not matrix:
            matrix = matrix_from_chars(characterlist)
        weight = up(
                tree,
                characterlist,
                matrix,
                pattern
                )
        root_weights[key] = min(weight[tree.root.name].values())
    return sum(root_weights.values())


        
