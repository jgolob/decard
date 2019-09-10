#!/usr/bin/env python
"""
        Goal: Take the target input, and covert it to something that can be dragged through the R
        Phyloseq package.
        Input: Targets.csv
        
"""
import argparse
import pandas as pd
from sklearn import linear_model

import json
import re
from Bio import SeqIO
import os

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
        args_parser.add_argument('OTU_test_dir', help='INPUT directory with otu test files')
        args_parser.add_argument('output_prefix', help="OUTPUT prefix to use (including path) for the two tables")
       
        
        args = args_parser.parse_args()
        
        source_dir = args.OTU_test_dir
        
        # Load data
        target_df = pd.read_csv(args.Target)
        
        
        # Convert weights into counts
        target_df['count'] = target_df.weight / target_df.weight.min()
        
        # Get our taxonomies
        target_df['taxonomy'] = target_df.apply(get_taxonomy_from_source,axis=1)
        target_df['taxonomy_list'] = target_df.taxonomy.apply(lambda t: t.split(';'))
        target_df['Order'] = target_df.taxonomy_list.apply(lambda tl: lookup_rank(tl,3))
        orders = target_df.Order.unique()
        
        
        by_community_sums = target_df.groupby(['community_name','taxonomy']).sum()
        by_community_sums_order = target_df.groupby(['community_name','Order']).sum()
        
        unique_taxonomies = target_df.taxonomy.unique()
        community_counts = pd.DataFrame(columns=unique_taxonomies)

        community_counts['Community'] = target_df.community_name.unique()
        community_counts.set_index('Community', inplace=True)
        
        community_counts_order = pd.DataFrame(columns=target_df.Order.unique())
        community_counts_order['Community'] = target_df.community_name.unique()
        community_counts_order.set_index('Community', inplace=True)
        
        for comm_tax, row in by_community_sums.iterrows():
                community_counts.loc[comm_tax[0], comm_tax[1]]=row['count']
        community_counts.fillna(0.0, inplace=True)
        
        for comm_order, row in by_community_sums_order.iterrows():
                community_counts_order.loc[comm_order[0], comm_order[1]]=row['count']
        community_counts_order.fillna(0.0, inplace=True)


        otu_test_fns = os.listdir(source_dir)
        otu_test_fns.sort()
        otu_test_fns.reverse()
        otu_test_fns = [fn for fn in otu_test_fns if fn[-4:]=='.csv']
        
        otu_tests = []
        for otu_test_fn in otu_test_fns:
            otu_test = pd.read_csv(os.path.join(source_dir,otu_test_fn))
            otu_test['raw_name'] = otu_test_fn
            otu_test['percent_dropped'] = otu_test.apply(lambda r: float(r.Dropped) / float(r.Total)*100.0, axis = 1)
            otu_test['label'] = otu_test.raw_name.apply(lambda rn: rn[:-12])
            otu_tests.append(otu_test)

        cumulative_df = pd.DataFrame()
        for otu_test in otu_tests:
                cumulative_df = cumulative_df.append(otu_test)
                
        longitude_df = pd.merge(cumulative_df, community_counts, left_on="Community", right_index=True)
        
        long_order_df = pd.merge(cumulative_df, community_counts_order, left_on="Community", right_index=True)
        
        
        # Regression time!
        clf = linear_model.LinearRegression()
        clf.fit(np.matrix([long_order_df.Bacteroidales]).T,long_order_df.specificity)
        clf.score(np.matrix([long_order_df.Bacteroidales]).T,long_order_df.specificity)
        
        
        
if __name__ == '__main__':
    main()
