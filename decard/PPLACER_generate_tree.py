#!/usr/bin/env python

import argparse
import pandas as pd
import sys
import json
from Bio import Phylo
import StringIO





def main():
    args_parser = argparse.ArgumentParser()

    args_parser.add_argument('jplace', help='jplace file from pplacer')
    args_parser.add_argument('otu_table', help='OTU table file for this jplace')
    args_parser.add_argument('tree_fn', help='Filname (and path) where to place the tree')

    args = args_parser.parse_args()
    
    otu_table = pd.read_csv(args.otu_table)
    
    jplace = json.load(open(args.jplace))
    # Extract tree, convert to string and put it in a stringIO object
    tree_h = StringIO.StringIO(str(jplace['tree']))
    # Load the tree
    tree = Phylo.read(tree_h,'newick')
    
    edge_num_index = jplace['fields'].index('edge_num')
    post_prob_idx = jplace['fields'].index('post_prob')
    distal_length_idx = jplace['fields'].index('distal_length')
    
    otu_id_to_edge = {}
    otu_ids = list(otu_table.otu_id.unique())
    for otu_id in otu_ids:
        for p_nm in jplace['placements']:
            placement_seq_ids = [nm[0] for nm in p_nm['nm'] ]
            if otu_id in placement_seq_ids:
                # sort by post_prob
                p_nm['p'].sort(key=lambda c: c[post_prob_idx])
                # As much as it pains me, take the top hit as that's what phyloseq will want. This is not great. 
                otu_id_to_edge[otu_id] = p_nm['p'][0]
                break
        
    
    """
        Multiple OTUs can map to the same edge, with some sort of distance from that edge
        Instead of renaming nodes, we need to add children.
    """
    for otu_id in otu_id_to_edge:
        changed = False
        for node in tree.find_elements("{"+str(otu_id_to_edge[otu_id][edge_num_index])+"}"):
            node.split(n=1, branch_length=otu_id_to_edge[otu_id][distal_length_idx])
            new_child = node.find_clades(node.name+"0").next()
            new_child.name = otu_id
            changed = True
        if not changed:
            sys.stderr.write("Failed to change "+str(otu_id)+" "+ str(otu_id_to_edge[otu_id][edge_num_index]+"\n"))
    
    Phylo.write(tree, args.tree_fn ,'newick')
    

if __name__ == '__main__':
    main()