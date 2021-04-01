#! /usr/bin/env python3.6

"""
Draw the tree and map multiples.

For multiple mapping, requires that all fragments of multiple have names ending
with '_1', '_2' and so on. If any leaf ID has this suffix, but doesn't match
other multiple, it is trimmed before processing further.
"""
import re

from argparse import ArgumentParser
from ete3 import Tree, TreeStyle, NodeStyle
from processing import change_support_format, trim_name, add_multi_annotation, \
    match_score, hmmer_name_mapping, are_ancestors

parser = ArgumentParser('Draw a tree and color multiples')
parser.add_argument('-t', type=str, help='Tree file (Newick)')
parser.add_argument('-s', type=float, default=0.5,
                    help='Match score threshold')
parser.add_argument('-o', type=str,
                    help='Image file. If not set, the tree is shown interactively.')
parser.add_argument('--support_threshold', default=0, type=float,
                    help='Collapse nodes below this support value')
parser.add_argument('--bracketed_support', action='store_true',
                    help='Assume support values are bracketed')
parser.add_argument('--quoted_names', action='store_true',
                    help='Assume leaf names are quoted')
parser.add_argument('--hmmer_ids', action='store_true',
                    help='Sequence IDs are produced by HMMER')
parser.add_argument('--skip_pairing', action='store_true',
                    help='Do not attempt to pair ancestral nodes, only mark leaves')
parser.add_argument('--circular', action='store_true',
                    help='Use circular tree style')
args = parser.parse_args()

################################################################################
# Preprocessing the tree line and the tree
################################################################################

with open(args.t) as tree_file:
    tree_line = tree_file.readline() # Assume there is a single line in the file
if args.bracketed_support:
    tree_line = change_support_format(tree_line)

# Preparing the tree for further analysis
tree = Tree(tree_line, quoted_node_names=args.quoted_names)
tree.set_outgroup(tree.get_midpoint_outgroup())

# Process raw names if necessary
if args.hmmer_ids:
    print('Processing HMMER IDs...')
    old_names = [x.name for x in tree.get_leaves()]
    name_map = hmmer_name_mapping(old_names)
    for leaf in tree.get_leaves():
        leaf.name = name_map[leaf.name]

for node in tree.traverse():
    if not node.is_leaf() and node.support < args.support_threshold:
        node.delete(prevent_nondicotomic=False, preserve_branch_length=True)

################################################################################
# Defining the multiple set
################################################################################

print('Selecting multiples...', end='')
multies = {} #Name-to-node mapping
multi_re = re.compile('(.+)_\d+$')

for leaf in tree.get_leaves():
    leaf.name = trim_name(leaf.name)
    if multi_re.match(leaf.name):
        multies[leaf.name] = leaf
print(f'{len(multies)} found')
# Trim subdomain numbers from non-multiple sequences
print('Trimming postfixes from non-multiples...')
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
# The dark gray bgcolor is set for all multiples. If a leaf is later assigned to
# the group, its style will be overridden. If not, this is left as the default.
multinode_style = NodeStyle()
multinode_style['bgcolor'] = 'gray'
for leaf in tree.get_leaves():
    if leaf in multies.values():
        leaf.set_style(multinode_style)

################################################################################
# Matching parts of a multidomain gene to each other and merging matches into
# match groups.
################################################################################

if args.skip_pairing:
    print('Not attempting to match ancestral nodes.')
else:
    print('Propagating multiples annotation...')
    # Mapping out descendants
    for node in tree.traverse(strategy='postorder'):
        add_multi_annotation(node, multies)
    # Selecting reciprocal best hit for each node with multidomain descendants

    # Do not process leaves to avoid bloating match set and performing costly
    # operations on them
    print('Generating match matrix...')
    node_refs = [x for x in tree.traverse('postorder') if x.multi_descendants and not x.is_leaf()]
    match_matrix = [[match_score(x.multi_descendants, y.multi_descendants) for x in node_refs]
                    for y in node_refs]
    matches = []
    print('Searching through the match matrix...', end='')
    for index, line in enumerate(match_matrix):
        for element in range(index+1, len(line)):
            if line[element] < args.s or match_matrix[element][index] < args.s:
                # Match quality cutoff. Can be tweaked to choose between larger
                # looser groups (at lower values) and exact matches (at 1.0)
                continue
            if (element, index) in matches:
                continue
            if are_ancestors(node_refs[index], node_refs[element]):
                continue
            matches.append((index, element))
    node_matches = [(node_refs[x[0]], node_refs[x[1]]) for x in matches]
    print(f'{len(node_matches)} raw matches found')
    print('Cleaning match set...', end='')
    # Preserve only the most ancient possible matches, discarding newer ones
    # In other words, if two clades match, remove all matches for their subclades
    clean_matches = []
    for match in node_matches:
        anc0 = match[0].get_ancestors()
        anc1 = match[1].get_ancestors()
        found = False
        for other_match in node_matches:
            if match == other_match:
                continue
            # Remove match (A,B) if there is match (C, D) such that C is the
            # ancestor of A and D is the ancestor of B
            if other_match[0] in anc0 and other_match[1] in anc1:
                found = True
                break
            elif other_match[1] in anc0 and other_match[0] in anc1:
                found = True
                break
            # Remove match (A, B) if there is match (A, C) such that C is the
            # ancestor of B
            if match[0] in other_match and (other_match[0] in anc1 or other_match[1] in anc1):
                found = True
                break
            if match[1] in other_match and (other_match[0] in anc0 or other_match[1] in anc0):
                found = True
                break
        if not found:
            clean_matches.append(match)

    tmp = set()
    for match in clean_matches:
        tmp.add(match[0])
        tmp.add(match[1])
    print(f'{len(clean_matches)} matches between {len(tmp)} nodes remaining')

    ############################################################################
    # Merge matches into match groups
    ############################################################################
    # For example, if three nodes have mutually overlapping sets of descendants,
    # they should become a single group of three, not three pairwise matches
    print('Merging matches into groups...', end='')
    groups = []
    for match in clean_matches:
        added = False
        for group in groups:
            if match[0] in group or match[1] in group:
                group.add(match[0])
                group.add(match[1])
                added = True
                break
        if not added:
            groups.append(set(match))

    print(f'produced {len(groups)} groups')
    # Perceptually distinct palette for groups. Copied from matplolib.cm.tab20
    colors = ('#1f77b4', '#aec7e8', '#ff7f0e', '#ffbb78', '#2ca02c',
              '#98df8a', '#d62728', '#ff9896', '#9467bd', '#c5b0d5',
              '#8c564b', '#c49c94', '#e377c2', '#f7b6d2', '#7f7f7f',
              '#c7c7c7', '#bcbd22', '#dbdb8d', '#17becf', '#9edae5')
    for i, group in enumerate(groups):
        # TODO: set leaf colors according to their groups
        tmp_style = NodeStyle()
        tmp_style['shape'] = 'circle'
        tmp_style['size'] = 15
        tmp_style['fgcolor'] = colors[i]
        tmp_style['bgcolor'] = colors[i]
        for node in group:
            # node.set_style(tmp_style)
            for leaf in node.get_leaves():
                leaf.set_style(tmp_style)

################################################################################
# Finally drawing the tree
################################################################################

tree_style = TreeStyle()
if args.circular:
    tree_style.mode = 'c'
tree_style.show_leaf_name = False
tree_style.show_branch_support = True
if args.o:
    tree.render(args.o, tree_style=tree_style)
else:
    tree.show(tree_style=tree_style)
