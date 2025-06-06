"""
Concordance factor computation.

This module provides basic code for the computation of the rooted site
concordance factor, an extension to the site concordance factor proposed by
Minh et al (2020). 

Rooted Site Concordance Factor
------------------------------

We introduce two major improvements over the site concordance factor discussed
by Minh et al. (2020). First, we provide a way to compute the site concordance
factor for rooted trees, introducing the notion of weakly and strongly
concordant quartets. Second, we apply significance testing in order to make
sure that the concordances obtained for a given tree and a given number of
characters are large enough to make it unlikely that they could have arisen by
chance. 

The distinction between weakly and strongly concordant quartets is derived from
the need that in order to be concordant with a rooted tree, not all quartets
need to be strict bipartitions that split the quartets into two. Alternatively,
they could also be what we call weakly concordant, and only separate the main
node being investigated during each step. 
"""

import statistics
import random
from tqdm import tqdm as progressbar
from pylostatistics.correlations import spearmanr, pointbiserialr
from pyloparsimony.util import pointbiserialr


def rooted_partition(tree, node):
    """
    Partition a tree into parts for the purpose of sampling.
    """
    # get the node and its descendants
    partA = [descendant.get_leaf_names() for descendant in
            tree[node].descendants]
    # get the next ancestral node
    siblings = [n for n in tree[node].ancestor.descendants if n.name !=
            node]
    # check if ancestral node is the root
    partB = [sibling.get_leaf_names() for sibling in siblings]
    if tree[node].ancestor.name == tree.root.name:
        if len(partB) > 1:
            return partA, partB
        return partA, [desc.get_leaf_names() for desc in
                siblings[0].descendants]
    else:
        visited_leaves = []
        for part in partA+partB:
            visited_leaves += part
        partB += [[name for name in tree.root.get_leaf_names() if name not in
            visited_leaves]]
        return partA, partB


def is_decisive(partitionA, partitionB):
    """
    Determine if a given quartet or other partition is decisive with respect to parsimony critiria.
    """
    characters = partitionA+partitionB
    counts = [characters.count(char) for char in set(characters)]
    if len(partitionA) in counts:
        return True
    return False


def is_concordant(partitionA, partitionB):
    """
    Determine concordance.

    .. note::
       
       Note that this returns concordance for structures [1123] and [1122], so
       the main criterion is that the first partition has the same character
       state which should not recur in the second partition.
    """
    if len(set(partitionA)) == 1 and partitionA[0] not in partitionB:
        return 1
    return 0


def select_leaves(partition, max_leaves=2):
    """
    Given a partition, select leaves which belong to this partition randomly.
    """
    out = []
    for part in partition:
        out += [random.choice(part)]
    return out[:max_leaves]


def select_chars(chars):
    """
    Select one character from a list of characters.
    """
    if isinstance(chars, (list, tuple)):
        return random.choice(chars)
    return chars


def randomly_concordant(chars):
    """
    This function estimates the concordance of a given quartet by chance.

    .. note::
       
       The estimate is based on the enumeration of all possible permutations of
       a decisive quartet. There are two basic structures, `1123` and `1122`.
       For the first structure, there are 12 permutations, with two permutations
       being concordant, and as a result, one out of six permutations will be
       concordant. For the second case, there are six possibilities, two of
       them concordant, and as a result, one out of three permutations is
       concordant.
       
    """
    assert len(chars) == 4
    if len(set(chars)) == 3:
        return 1/6
    elif len(set(chars)) == 2:
        return 1/3


#def is_random_concordant(charsA, charsB):
#    """
#    Function checks if a shuffled version of the partitions is concordant.
#    """
#    chars = charsA + charsB
#    random.shuffle(chars)
#    charsAn, charsBn = chars[:len(charsA)], chars[len(charsA):]
#    return is_concordant(charsAn, charsBn)



def rooted_site_concordance(
        tree, 
        patterns, 
        iterate=20, 
        missing="Ø",
        max_leaves=2, 
        iterate_correlation=100, 
        correlation=pointbiserialr,
        ):
    """
    Rooted site concordance factor.

    .. note::
       
       The output is a dictionary in which the node labels from the tree are
       the keys, and each node has additional information for each pattern in
       the original data, given in the form of a list:
       
       * "attested" refers to the sites and provides the mean score of all trials
         the site concordance
       * "expected" refers to the expected concordance values of the sites
       * "decisive" indicates for each site if it is decisive or not
       * "patterns" provides the pattern identifiers ("cognate IDs")
       * "trials" provides the information on the number of effective trials
         which were carried out for each site
       * chars provides detailed information for all characters
       
       Concordance factors are available in two more keys:

       * "concordance" gives the concordance factors for all
         sites that were decisive
       * "concordance_expected" gives the expected average concordance for a
         node
       * "r" gives information on the correlation
       * "p" is the p-value of the correlation
       * "oddsratio" is the odds ratio computed from the attested and the
         expected concordance scores. 
    """
    nodes = {
            node.name: {
                "attested": [], 
                "chars": [], 
                "expected": [], 
                "decisive": [], 
                "patterns": [], 
                "trials": []
                } for node in tree.preorder[1:] if node.descendants}
    for node in progressbar(tree.preorder[1:], desc="computing concordance"):
        if node.descendants:
            partA, partB = rooted_partition(tree, node.name)
            for pid, pattern in patterns.items():
                attested, expected, visited, chars = [], [], set(), []
                for i in range(iterate):
                    choiceA, choiceB = (
                            select_leaves(partA, max_leaves=max_leaves),
                            select_leaves(partB, max_leaves=max_leaves)
                            )
                    if tuple(choiceA+choiceB) not in visited and len(choiceA+choiceB) >= 4:
                        visited.add(tuple(choiceA+choiceB))
                        charsA, charsB = (
                                [select_chars(pattern[taxon]) for taxon in choiceA], 
                                [select_chars(pattern[taxon]) for taxon in choiceB]
                                )
                        if not missing in charsA+charsB:
                            if is_decisive(charsA, charsB):
                                a = is_concordant(charsA, charsB)
                                attested += [a]
                                expected += [randomly_concordant(charsA + charsB)]
                                chars += [(charsA, charsB, a,
                                    expected[-1])]
                if attested:
                    nodes[node.name]["chars"] += [chars]
                    nodes[node.name]["attested"] += [sum(attested)]
                    nodes[node.name]["decisive"] += [len(attested)]
                    nodes[node.name]["patterns"] += [pid]
                    nodes[node.name]["trials"] += [len(visited)]
                    nodes[node.name]["expected"] += [sum(expected)]
    for node in progressbar(nodes, desc="checking significance"):
        decisive_a = [n / d for n, d in zip(nodes[node]["attested"], nodes[node]["decisive"]) if d]
        decisive_e = [n / d for n, d in zip(nodes[node]["expected"], nodes[node]["decisive"]) if d]
        if decisive_a:
            sites = len(nodes[node]["attested"])
            
            # compute the correlation to test the significance
            lstA = decisive_a + decisive_e
            lstB = [1 for x in decisive_a] + [0 for x in decisive_e]
            try:
                r, p = correlation(lstA, lstB, iterate=iterate_correlation)
            except ZeroDivisionError:
                r, p = 0, 1
                print(lstA, lstB)
            c_len = len([x for x, y in zip(decisive_a, decisive_e) if x >= y])
            d_len = len([x for x, y in zip(decisive_a, decisive_e) if x < y])
            nodes[node]["oddsratio"] = \
                    statistics.mean(decisive_a) / statistics.mean(decisive_e)
            nodes[node]["concordance"] = statistics.mean(decisive_a)
            nodes[node]["concordance_expected"] = statistics.mean(decisive_e)
            nodes[node]["r"] = r
            nodes[node]["p"] = p
        else:
            nodes[node]["r"] = 0
            nodes[node]["oddsratio"] = 0
            nodes[node]["p"] = 1
            nodes[node]["concordance"] = 0
            nodes[node]["concordance_expected"] = 0
    return nodes



def concordance_statistics(nodes, significant=0.05):
    out = {"significant": 0, "nodes": len(nodes), "significant_nodes": []}
    for node, items in nodes.items():
        if items["p"] <= significant:
            out["significant"] += 1
            out["significant_nodes"] += [(items["oddsratio"], items["r"])]
    out["significant_proportion"] = out["significant"] / out["nodes"]
    if out["significant_nodes"]:
        out["significant_oddsratio"] = statistics.mean(
                [x[0] for x in out["significant_nodes"]])
        out["significant_correlation"] = statistics.mean(
                [x[1] for x in out["significant_nodes"]])
    else:
        out["significant_nodes"] = 0
        out["significant_oddsratio"] = 0
        out["significant_r"] = 0
    return out


