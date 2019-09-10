#!/usr/bin/env python
"""
    Take a biom format file and flattens it into a table format needed for subsequent analysis
    The existing utility is very memory intensive. This is meant to be a lightweight implementation

"""

import argparse
import pandas as pd
import biom
import csv
import sys



def sparse_matrix_to_otu_index(sm):
    """
        Raw data in many of these formats is sparse, using numpy's sparse matrix
        Within the spase matrix in the [1][0]
        
        Input: A numpy sparce matrix
        Output: list of OTU matrix Ids and weights
        
    """
    try:
        return sm.nonzero()[1][0]
    except IndexError:
        return None

def get_weight(row):
    try:
        return row.raw[0,row.otu_index]
    except ValueError:
        return 0
    
def main():
    args_parser = argparse.ArgumentParser()

    args_parser.add_argument('input', help="biom (v2) format input")
    args_parser.add_argument('output', help="TSV file into which to place our table")

    args = args_parser.parse_args()
    
    out_h = open(args.output,'w')
    
    writer = csv.writer(out_h, delimiter='\t')    
    
    # Load the biom table and extract the seq IDs and otu IDs
    biom_table = biom.load_table(args.input)
    seq_ids = biom_table.ids('sample')
    # Seq ids need to be cleaned up
    seq_ids = [seq_id.replace('#OTU ID',"") for seq_id in seq_ids]
    
    otu_ids = biom_table.ids('observation')
    
    
    # Make some headers for our table. First line is rote text
    writer.writerow(['# Constructed from biom file\n'])
    
    # second line is all the seq ids, with a rote first column
    writer.writerow(['#OTU ID']+seq_ids)
    
    # Use list comprehension to extract the otu_index and weight for each sequence
    # The index of this list will be the same as the seq_ids list above
    # Still relatively lightweight, as it's just a list of tuples of an int and float
    # This is the longest and most memory intensive task of the process.
    seq_to_otu = [(
        raw.nonzero()[1][0], # otu_index
        raw[0,raw.nonzero()[1][0]] # weight
    )
    for raw in biom_table.iter_data(dense=False)]
    
    # Convert to a pandas DF to help us transpose to have OTU_ids be the row, and output our data
    seq_to_otu_df = pd.DataFrame(seq_to_otu, columns=['otu_index','weight'])
    for otu_index, otu_block in seq_to_otu_df.groupby('otu_index'):
        sys.stdout.write("OTU: "+ str(otu_ids[otu_index])+"\n")
        # The row is the otu ID followed by the weight for each seq, or zero if there is no entry.
        # We'll use list comprehension to do this somewhat quickly
        writer.writerow([otu_ids[otu_index]]+[otu_block.weight.loc[seq_index]  if seq_index in otu_block.index else 0.0 for seq_index in xrange(len(seq_ids))])
    
    out_h.close()
    
if __name__ == '__main__':
    main()