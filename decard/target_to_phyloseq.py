#!/usr/bin/env python
"""
        Goal: Take the target input, and covert it to something that can be dragged through the R
        Phyloseq package.
        Input: Targets.csv
        
"""
import argparse
import pandas as pd

import json
import re
from Bio import SeqIO

header_re = re.compile('description=\"(?P<description>[^\"]+)\" organism=\"(?P<organism>[^\"]+)\" taxonomy=\"(?P<taxonomy>[^\"]+)\" ncbi_tax_id=\"(?P<ncbi_tax_id>[^\"]+)\"')

def get_taxonomy_from_source(row):
        try:
                seq_recs = SeqIO.parse(row.source_file,'fasta')
        
                for sr in seq_recs:
                    if sr.id == row.sequence_id:
                        m = header_re.search(sr.description)
                        if m:
                            return m.group('taxonomy')
                        else:
                            return ""
                        break
                if not 'ncbi_tax_id' in target:
                    return None
        except IOError:
                return ""
def lookup_rank(taxonomy_list, rank_index):
        try:
                return taxonomy_list[rank_index]
        except IndexError:
                return None

def main():
        args_parser = argparse.ArgumentParser()
        
        args_parser.add_argument('Target', help='INPUT  (target) file, in CSV format')
        args_parser.add_argument('output_prefix', help="OUTPUT prefix to use (including path) for the two tables")
       
        
        args = args_parser.parse_args()
        
        
        # Load data
        target_df = pd.read_csv(args.Target)
        
        
        # Convert weights into counts
        target_df['count'] = target_df.weight / target_df.weight.min()
        
        # Get our taxonomies
        target_df['taxonomy'] = target_df.apply(get_taxonomy_from_source,axis=1)
        target_df['taxonomy_list'] = target_df.taxonomy.apply(lambda t: t.split(';'))
        
        # Start in on our matrix
        
        # Get our source IDS
        source_IDs = target_df.sequence_id.unique()
        
        otu_matrix = pd.DataFrame(index=source_IDs)
        
        for comm_name, comm_block in target_df.groupby('community_name'):
                otu_matrix[comm_name] = [round(comm_block[comm_block.sequence_id == otu_id]['count'].sum()) for otu_id in otu_matrix.index]
        
        # Now our tax table
        
        tax_table = pd.DataFrame(index=source_IDs, columns=['Domain', 'Phylum', 'Class', 'Order','Family','Genus', 'Species'])

        # This isn't very pandas, but it'll work for now
        for sequence_id, row in tax_table.iterrows():
                #print sequence_id, len(target_df[target_df.sequence_id==sequence_id])
                taxonomy_list = target_df[target_df.sequence_id==sequence_id].iloc[0].taxonomy_list
                species = target_df[target_df.sequence_id==sequence_id].iloc[0].species
                # Populate
                row.Domain = lookup_rank(taxonomy_list,0)
                row.Phylum = lookup_rank(taxonomy_list,1)
                row.Class = lookup_rank(taxonomy_list,2)
                row.Order = lookup_rank(taxonomy_list,3)
                row.Family = lookup_rank(taxonomy_list,4)
                row.Genus = lookup_rank(taxonomy_list,5)
                row.Species = species
                
        
        # OUTPUT TIME!
        tax_table.to_csv(args.output_prefix+'.tax_table.csv')
        otu_matrix.to_csv(args.output_prefix+".otu_matrix.csv")

if __name__ == '__main__':
    main()
