#!/usr/bin/python

import sys

input=sys.argv[1]
output=sys.argv[2]
command=sys.argv[3]

tree_file=open(input,"r")
tree_output=open(output,"w")

#lines=tree_file.readlines()

new_line = list()
for line in tree_file:
	if line.startswith("\tTREE"):
		new_line.append(line)
	else:
		pass
	
out_line = list()
tree_name_line = list()

tree_output.write("#NEXUS\n")
tree_output.write("BEGIN TREES;\n")

for i in xrange(len(new_line)):
	TreeName="tree_{0}".format(i+1)
	h=new_line[i].split(" ")
	out_line=[h[0],TreeName,h[3],h[4],h[5]]
	tree_name_line.append(TreeName)
	tree_output.write(' '.join(out_line))
	

tree_output.write("END;")
tree_output.write("\n")
tree_output.write("BEGIN PHYLONET;\n")
tree_output.write("\t")
tree_output.write(''.join(command))
tree_output.write(";\n")
tree_output.write("END;")
tree_output.write("\n")
tree_output.close()
