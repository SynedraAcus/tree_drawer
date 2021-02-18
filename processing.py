"""
Various functions for tree and name preprocessing
"""

import re


def sub_replacement(match):
    return match.group(2)+':'+match.group(1)


def change_support_format(tree_line):
    """
    Change support value format from bracketed to ete3-compatible colons

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
    :param name:
    :return:
    """
    # TODO: Maybe there is something else that needs trimming?
    return re.sub('_\(\d+-\d+\)', lambda x: '', name)


def add_multi_annotation(node, multies):
    """
    Add a list of multiples descending from this node to node annotation
    :param node:
    :return:
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
    Match score between prefix lists
    :param vector1:
    :param vector2:
    :return:
    """
    # Otherwise, a raw overlap
    s1 = set(vector1)
    s2 = set(vector2)
    # if len(s1) < len(vector1) or len(s2) < len(vector2):
    #     # If a node is a descendant to both parts of some sequence, it is not
    #     # what we're interested in
    #     #
    #     # This part may screw us over if there is a lot of nested duplications
    #     return 0
    # else:
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
