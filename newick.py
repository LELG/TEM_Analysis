from ete3 import Tree, TreeStyle, TextFace
#t = Tree( "((5),0,2,3);", format=1 )
#t = Tree( "(3,5);", format =1 )
t = Tree("(A:1,(B:1,(E:1,D:1):0.5):0.5);" )
t = Tree( "( ( A,B,(5) ) );" ,format=1)
t = Tree.from_parent_child_table([("A", "B"), ("A", "C"), ("C", "D"), ("C", "E")])
t = Tree.from_parent_child_table([('P', 'PC1'), ('PC1', 'PC1,2-26:4402'), ('PC1', 'PC1,3'), ('PC1,3', 'PC1,3,5-26:6743'), ('PC1', 'PC1,14'), ('PC1,3,14', 'PC1,3,14,4-27:124'), ('PC1', 'PC1,9'), ('PC1,9', 'PC1,9,7-27:675'), ('PC1', 'PC1,0-26:4132')])

t = Tree.from_parent_child_table([("P", "PC0-23:7392", 100), ("P", "PC1-25:6592", 7.5), ("PC1-25:6592", "PC1,1-26:7008", 6), ("PC1-25:6592", "PC1,2-26:7008", 6)])
ts = TreeStyle()
ts.show_leaf_name = True
#ts.title.add_face(TextFace("Hello ETE", fsize=20), column=0)
t.show(tree_style=ts)