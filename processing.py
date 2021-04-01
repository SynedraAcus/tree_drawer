"""
Various functions for tree and name preprocessing
"""

import re


def sub_replacement(match):
    return match.group(2)+':'+match.group(1)


def change_support_format(tree_line):
    """
    Change Newick branch support format from Biopython-compatible
    bracketed to ete3-compatible colons

    For example:
    ((A:1, B:0.7)0.8[65], C) becomes ((A:1, B:0.7)65:0.8, C)
    :param str tree_line: Parseable Newick
    :return:
    """
    node_re = re.compile(':([\d.E-]+)\[([\d.E-]+)\]')
    return re.sub(node_re, sub_replacement, tree_line, count=0)


def trim_name(name):
    """
    Trim the leaf name.

    Remove domain position (if any)
    """
    return re.sub('_\(\d+-\d+\)', lambda x: '', name)


def add_multi_annotation(node, multies):
    """
    Add a list of multiples descending from this node to node annotation

    Assumes leaf-first traversal. If called for a node before any of its
    descendants, this function will break.
    """
    if node.is_leaf():
        # For a leaf, add either name or nothing
        node.add_feature('multi_descendants', [])
        if node.name in multies:
            node.multi_descendants.append('_'.join(node.name.split('_')[:-1]))
    else:
        node.add_feature('multi_descendants', [])
        for child in node.children:
            node.multi_descendants += child.multi_descendants


def match_score(vector1, vector2):
    """
    Match score between prefix lists.
    :param vector1:
    :param vector2:
    :return:
    """
    s1 = set(vector1)
    s2 = set(vector2)
    return len(s1.intersection(s2))/len(s1.union(s2))


def maxindices(l):
    """
    Get indices for all occurences of maximal element in list
    :param l:
    :return:
    """
    max_indices = []
    max_value = l[0] #Assume un-exhaustible iterator
    for i, v in enumerate(l):
        if v > max_value:
            max_value = v
            max_indices = [i]
        elif v == max_value:
            max_indices.append(i)
    return max_indices


def hmmer_name_mapping(names):
    """
    Turn HMMER-produced list of IDs into numbered and trimmed leaf IDs
    Takes a list of names, returns a dict {old_name: new_name}
    """
    r = {}
    for name in names:
        # Remove second part of name, it's useless
        r[name] = name.split(' [subseq')[0]
    used_queries = set()
    first_pos_re = re.compile('/(\d+)-\d+ ')
    for name in r:
        # Renumber subdomains so that they use their position
        # in a given protein instead of coordinates
        if r[name][-2] == '_':
            # If the name has already been postfixed, it cannot be processed
            # by the code below
            continue
        query = r[name].split('/')[0]
        if query in used_queries:
            continue
        used_queries.add(query)
        matches = {x: r[x] for x in r if query in x}
        if len(matches) > 1:
            # Add correct number
            positions = {}
            for x in matches:
                re_match = first_pos_re.search(x)
                first_pos = re_match.groups(1)
                positions[x] = int(first_pos[0])
            first_pos = list(positions.values())
            first_pos.sort()
            for name in matches:
                r[name] = query + '_' + str(first_pos.index(positions[name]) + 1)
        else:
            # Don't add postfixes and discard coordinates
            # if it is the only domain in a sequence
            r[name] = query + '_1'
    return r


def are_ancestors(node1, node2):
    """
    Take two ete3 tree nodes.

    Return True if either is the ancestor of another
    :param node1:
    :param node2:
    :return:
    """
    for ancestor in node1.iter_ancestors():
        if ancestor == node2:
            return True
    for ancestor in node2.iter_ancestors():
        if ancestor == node1:
            return True
    return False