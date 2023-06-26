from ete3 import Tree
import os
import sys

def get_dist(tree:str):  
    #'SpeciesTree_rooted.txt'
    ORTHO_TREE = os.path.join(tree)
    
    dict_tree = dict()
    tree = Tree(ORTHO_TREE)
    for sp1 in tree.iter_leaves():
        dict_dist = dict()
        for sp2 in tree.iter_leaves():
            if sp1 != sp2:
                dict_dist[str(sp2)[3:]] = sp1.get_distance(sp2)
        dict_tree[str(sp1)[3:]] = dict_dist

    return dict_tree
    

if __name__ == "__main__":
    get_dist(sys.argv[1])
