#from ete3 import Tree
import sys

# parse tree
with open(sys.argv[1], 'r') as nwk:
	nwk_string = ''
	for line in nwk:
		nwk_string += line.rstrip("\n")

tree = Tree(nwk_string)

# label nodes
node_num = 0
for node in tree.traverse("postorder"):
	if len(node.name) == 0:
		node.add_features(name=str(node_num))
		node_num += 1

#print tree.get_ascii(attributes=["name"], show_internal=True)

# parse list of leaves of interest
with open(sys.argv[2], 'r') as list:
	leaf_list = []
	for line in list:
		leaf_list.append(line.rstrip("\n"))

# node of interest
interest_node_name = sys.argv[3]
for node in tree.traverse("postorder"):
	if node.name == interest_node_name:
		interest_node = node

# print branch lengths
for leaf in leaf_list:
	print leaf + "\t" + str(tree.get_distance(leaf, interest_node))
