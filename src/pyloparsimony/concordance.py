"""
Concordance factor computation.
"""
import statistics
import random
from tqdm import tqdm as progressbar

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
    characters = partitionA+partitionB
    counts = [characters.count(char) for char in set(characters)]
    if len(partitionA) in counts:
        return True
    return False


def is_concordant(partitionA, partitionB):
    if len(set(partitionA)) == 1 and partitionA[0] not in partitionB and len(set(partitionB)) == 1:
        return 1
    return 0


def select_leaves(partition):
    out = []
    for part in partition:
        out += [random.choice(part)]
    return out


def select_chars(chars):
    if isinstance(chars, (list, tuple)):
        return random.choice(chars)
    return chars


def rooted_site_concordance(tree, patterns, iterate=20, missing="Ø"):
    """
    Rooted site concordance factor.

    .. note::
       
       The output is a dictionary in which the node labels from the tree are
       the keys, and each node has additional information for each pattern in
       the original data, given in the form of a list:
       
       * "sites" refers to the sites and provides the mean score of all trials
         the site concordance
       * "decisive" indicates for each site if it is decisive or not
       * "patterns" provides the pattern identifiers ("cognate IDs")
       * "trials" provides the information on the number of effective trials
         which were carried out for each site
       
       Concordance factors are available in two more keys:

       * "concordance_all_sites" gives information for the concordance factor
         computed from all sites, even a site yielded only trials that were not
         decisive
       * "concordance_decisive_sites" gives the concordance factors for all
         sites that were decisive
    """
    nodes = {node.name: {"sites": [], "rscores": [], "decisive": [], "patterns": [], "trials": []} for node in tree.preorder[1:] if
            node.descendants}
    for node in progressbar(tree.preorder[1:], desc="computing concordance"):
        if node.descendants:
            partA, partB = rooted_partition(tree, node.name)
            for pid, pattern in patterns.items():
                score, rscore, visited = [], [], set()
                for i in range(iterate):
                    choiceA, choiceB = (
                            select_leaves(partA),
                            select_leaves(partB)
                            )
                    if tuple(choiceA+choiceB) not in visited and len(choiceA+choiceB) >= 4:
                        visited.add(tuple(choiceA+choiceB))
                        charsA, charsB = (
                                [select_chars(pattern[taxon]) for taxon in choiceA], 
                                [select_chars(pattern[taxon]) for taxon in choiceB]
                                )
                        if not missing in charsA+charsB:
                            if is_decisive(charsA, charsB):
                                score += [is_concordant(charsA, charsB)]
                                rscore += [is_random_concordant(charsA,
                                    charsB)]
                if score:
                    nodes[node.name]["sites"] += [statistics.mean(score)]
                    nodes[node.name]["decisive"] += [len(score)]
                    nodes[node.name]["patterns"] += [pid]
                    nodes[node.name]["trials"] += [len(visited)]
                    nodes[node.name]["rscores"] += [statistics.mean(rscore)]
                else:
                    nodes[node.name]["sites"] += [0]
                    nodes[node.name]["decisive"] += [0]
                    nodes[node.name]["patterns"] += [pid]
                    nodes[node.name]["trials"] += [len(visited)]
                    nodes[node.name]["rscores"] += [0]
    for node in nodes:
        nodes[node]["concordance_all_sites"] = statistics.mean(nodes[node]["sites"])
        nodes[node]["random_concordance_all_sites"] = statistics.mean(
                nodes[node]["rscores"])
        decisive = [n for n, d in zip(nodes[node]["sites"], nodes[node]["decisive"]) if d]
        rdecisive = [n for n, d in zip(nodes[node]["rscores"], nodes[node]["decisive"]) if d]
        if decisive:
            nodes[node]["concordance_decisive_sites"] = statistics.mean(decisive)
            nodes[node]["random_concordance_decisive_sites"] = statistics.mean(rdecisive)
            try:
                nodes[node]["corrected_concordance"] = \
                        nodes[node]["concordance_decisive_sites"] / nodes[node][
                                "random_concordance_decisive_sites"]
            except ZeroDivisionError:
                nodes[node]["corrected_concordance"] = 0

        else:
            nodes[node]["concordance_decisive_sites"] = 0
            nodes[node]["random_concordance_decisive_sites"] = 0
            nodes[node]["corrected_concordance"] = 0
    return nodes


def is_random_concordant(charsA, charsB):
    chars = charsA + charsB
    random.shuffle(chars)
    charsAn, charsBn = chars[:len(charsA)], chars[len(charsA):]
    return is_concordant(charsAn, charsBn)



