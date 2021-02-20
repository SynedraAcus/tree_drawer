#! /usr/bin/env python3

"""
Draw the tree and map multiples.

For multiple mapping, requires that all fragments of multiple have names ending
with '_1', '_2' and so on. If any leaf ID has this suffix, but doesn't match
other multiple, it is trimmed before processing further.
"""
import re

from argparse import ArgumentParser
from collections import defaultdict
from functools import reduce
from ete3 import Tree, TreeStyle, NodeStyle
from processing import change_support_format, trim_name, add_multi_annotation, \
    match_score, maxindices

parser = ArgumentParser('Draw a tree and color multiples')
parser.add_argument('-t', type=str, help='Tree file (Newick)')
parser.add_argument('-s', type=float, default=-0.4,
                    help='Match score threshold')
parser.add_argument('--bracketed_support', action='store_true',
                    help='Assume support values are bracketed')
parser.add_argument('--quoted_names', action='store_true',
                    help='Assume leaf names are quoted')
args = parser.parse_args()

# Preprocessing the tree line
with open(args.t) as tree_file:
    tree_line = tree_file.readline() # Assume there is a single line in the file
if args.bracketed_support:
    tree_line = change_support_format(tree_line)

tree = Tree(tree_line, quoted_node_names=args.quoted_names)
tree.set_outgroup(tree.get_midpoint_outgroup())
# Defining the multiple set
multies = {} #Name-to-node mapping
multi_re = re.compile('(.+)_\d+$')

for leaf in tree.get_leaves():
    leaf.name = trim_name(leaf.name)
    if multi_re.match(leaf.name):
        multies[leaf.name] = leaf
# Trim subdomain numbers from non-multiple sequences
# This part is ugly, but it works for now
to_trim = {} #Name:prefix
for name in multies:
    prefix = '_'.join(name.split('_')[:-1])
    matches = [x for x in multies if prefix in x]
    if len(matches) < 2:
        to_trim[name] = prefix
for name in to_trim:
    #multies[to_trim[name]] = multies[name]
    multies[name].name = to_trim[name]
    del multies[name]
# Set the markers right now, but leave colors for when we have clades
multinode_style = NodeStyle()
multinode_style['shape'] = 'circle'
multinode_style['size'] = 10
for leaf in multies.values():
    leaf.set_style(multinode_style)

# Mapping out descendants
for node in tree.traverse(strategy='postorder'):
    add_multi_annotation(node, multies)

################################################################################
# Selecting reciprocal best hit for each node with multidomain descendants
################################################################################

# Do not process leaves to avoid bloating match set and performing costly
# operations on them
node_refs = [x for x in tree.traverse() if x.multi_descendants and not x.is_leaf()]
match_matrix = [[match_score(x.multi_descendants, y.multi_descendants) for x in node_refs]
                for y in node_refs]
matches = []
for index, line in enumerate(match_matrix):
    for element in range(len(line)):
        # Could use less checks, but at least this is readable
        if line[element] < args.s:
            # Match quality cutoff. The exact value is arbitrary, but anything
            # above 0.01 produces the same result in Chalcone
            continue
        if element == index:
            # No self matches
            continue
        if match_matrix[element][index] < args.s:
            continue
        # if index not in maxindices(match_matrix[element]):
        #     # Not reciprocal
        #     continue
        if node_refs[index] in node_refs[element].iter_descendants():
            # Ancestor and descendant should not match each other
            continue
        if node_refs[element] in node_refs[index].iter_descendants():
            continue
        if ((element, index)) not in matches:
            # This match haven't been found before
            matches.append((index, element))
node_matches = [(node_refs[x[0]], node_refs[x[1]]) for x in matches]
clean_matches = []
while node_matches:
    match = node_matches.pop()
    anc1 = match[0].get_ancestors()
    anc2 = match[1].get_ancestors()
    found = False
    for other_match in node_matches:
        # Remove match (A,B) if there is match (C, D) such that C is the
        # ancestor of A and D is the ancestor of B
        if other_match[0] in anc1 and other_match[1] in anc2:
            found = True
            break
        elif other_match[1] in anc1 and other_match[0] in anc2:
            found = True
            break
        # Remove match (A, B) if there is match (A, C) such that C is the
        # ancestor of B
        if match[0] in other_match and (other_match[0] in anc2 or other_match[1] in anc2):
            found = True
            break
        if match[1] in other_match and (other_match[0] in anc2 or other_match[1] in anc2):
            found = True
            break
    if not found:
        clean_matches.append(match)

# Random debug shit
tmp = set()
for match in clean_matches:
    tmp.add(match[0])
    tmp.add(match[1])
print(f'Found {len(clean_matches)} matches between {len(tmp)} nodes')
# for match in clean_matches:
#     print('MATCH')
#     print(sorted([x.name for x in match[0].get_leaves()]))
#     print(sorted([x.name for x in match[1].get_leaves()]))

################################################################################
# Merge matches into match groups
################################################################################
# For example, if three nodes have mutually overlapping sets of descendants,
# they should become a single group of three, not three pairwise matches
#TODO: fix group generation
groups = []
for match in clean_matches:
    added = False
    for group in groups:
        if match[0] in group or match[1] in group:
            group.add(match[0])
            group.add(match[1])
            added = True
            continue
    if not added:
        groups.append(set(match))
# while clean_matches:
#     group = set(clean_matches.pop())
#     to_remove = []
#     for other_match in clean_matches:
#         if other_match[0] in group or other_match[1] in group:
#             group = group.union(set(other_match))
#     if group not in groups:
#         groups.append(group)
#     # by_node[match[0]] += match
#     # by_node[match[1]] += match
for group in groups:
    for other_group in groups:
        if group.intersection(other_group):
            print('OVERLAP')

print(f'They were merged into {len(groups)} groups')
colors = ('#1f77b4', '#aec7e8', '#ff7f0e', '#ffbb78', '#2ca02c',
          '#98df8a', '#d62728', '#ff9896', '#9467bd', '#c5b0d5',
          '#8c564b', '#c49c94', '#e377c2', '#f7b6d2', '#7f7f7f',
          '#c7c7c7', '#bcbd22', '#dbdb8d', '#17becf', '#9edae5')
for i, group in enumerate(groups):
    tmp_style = NodeStyle()
    tmp_style['shape'] = 'circle'
    tmp_style['size'] = 15
    tmp_style['fgcolor'] = colors[i]
    for node in group:
        node.set_style(tmp_style)

# Drawing the Tree
tree_style = TreeStyle()
# tree_style.mode = 'c'
tree.show(tree_style=tree_style)
