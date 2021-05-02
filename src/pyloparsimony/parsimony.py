"""
Compute parsimony analyses

"""

from pylotree import Tree
from collections import defaultdict
from itertools import product

def parsimony_up(
        patterns=None,
        taxa=None,
        tree=None,
        matrix=None,
        characters=None,
        missing="Ø",
        output=None):
    """
    Carries out sankoff parsimony.
    """
    W = defaultdict(dict)
    lookupT = {name: i for i, name in enumerate(taxa)}
    lookupM = {char: i for i, char in enumerate(characters)}

    if not isinstance(tree, Tree):
        tree = Tree(tree)

    for node, name in map(lambda x: (x, x.name), tree.postorder):        
        if not node.descendants:
            for i, char in enumerate(characters):
                if char in patterns[lookupT[name]]:
                    W[name][char] = 0
                else:
                    W[name][char] = 1000000
        else:
            for i, charA in enumerate(characters):
                scoresA = []
                for nodeB, nameB in map(
                        lambda x: (x, x.name),
                        tree[name].descendants):
                    scoresB = []
                    for j, charB in enumerate(characters):
                        weightB = W[nameB][charB]
                        weightNew = weightB+matrix[i][j]
                        scoresB += [weightNew]
                    smin = min(scoresB)
                    scoresA += [smin]
                W[name][charA] = sum(scoresA)
    
    if output == "weight":
        return min(W[tree.name].values())
    if output == "chars":
        return [x for x, y in W[tree.root].items() if y ==
                min(W[tree.name].values())]
    return W


def parsimony_down(
        weights=None,
        patterns=None,
        taxa=None,
        tree=None,
        matrix=None,
        characters=None
        ):
    if not isinstance(tree, Tree):
        tree = Tree(tree)
    smin = min(weights[tree.name].values())
    root_chars = [a for a, b in weights[tree.name].items() if b == smin]
    
    # prepare the queue
    queue = []
    for char in root_chars:
        nodes = []
        for child in tree[tree.name].descendants:
            nodes += [child]
        queue += [([(nodes, tree.name, char, characters.index(char))], [(tree.name, char)])]
    
    output = []
    while queue:
        nodes, scenario = queue.pop(0)
        if not nodes:
            output += [scenario]
        else:
            children, parent, pchar, pidx = nodes.pop()
            pscore = weights[parent][pchar]
            for comb in product(*len(children)*[characters]):
                score = 0
                for i, char in enumerate(comb):
                    cidx = characters.index(char)
                    score += matrix[pidx][cidx]
                    score += weights[children[i].name][char]
                if score == pscore:
                    new_nodes = [n for n in nodes]
                    new_scenario = [s for s in scenario]
                    for child, char in zip(children,comb):
                        new_nodes += [(tree[child.name].descendants, child.name, char,
                            characters.index(char))]
                        new_scenario += [(child.name, char)]
                    queue += [(new_nodes, new_scenario)]
    return output

def parsimony(
        patterns=None,
        taxa=None,
        tree=None,
        matrix=None,
        characters=None
        ):
    if not isinstance(tree, Tree):
        tree = Tree(tree)

    W = parsimony_up(
            patterns=patterns,
            taxa=taxa,
            tree=tree,
            matrix=matrix,
            characters=characters,
            )
    
    # get minimal weight
    smin = min(W[tree.name].values())
    weights = [b for a, b in W[tree.name].items() if b == smin]
    scenarios = parsimony_down(
            weights=W,
            patterns=patterns,
            taxa=taxa,
            tree=tree,
            matrix=matrix,
            characters=characters
            )
    
    return weights, scenarios, W



