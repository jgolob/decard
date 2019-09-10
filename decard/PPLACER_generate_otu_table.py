#!/usr/bin/env python

import argparse
import sqlite3
import pandas as pd
import sys




def main():
    args_parser = argparse.ArgumentParser()

    args_parser.add_argument('--dedup_info','-di', help='Headerless CSV file with a mapping of OTU to representitive sequence from each sample or a sequence', required=True)
    args_parser.add_argument('--seq_info', '-s', help='Headerless CSV file with a mapping of sequence to a sample')
    args_parser.add_argument('--placements_db','-p', help='SQLite3 database with the placement outputs', required=True)
    args_parser.add_argument('--bySeqOutput','-o', help="CSV file for output", required=True)
    args_parser.add_argument('--likelihood_threshold', '-lt', default = 0.95, help='Minimum likelihood to make a call')
    args_parser.add_argument('--clusters', '-c', help='OPTIONAL input. A whitespace delimited file where the first column is the OTU rep ID, the second column the dereplication ID and the third a common seperated list of all member sequences')

    args = args_parser.parse_args()
    
    dedup_info_df = pd.read_csv(args.dedup_info, header=-1, names=['otu_id','seq','weight'])
    seq_info_df = pd.read_csv(args.seq_info, header=-1, names=['seq','community'])
    
    if args.clusters:
        # Load the cluster info into a dataframe
        clusters_df = pd.read_csv(args.clusters, header=-1, names=['otu_id', 'derep_id', 'seq_ids_raw'], delimiter=" ")
        # Split out the seq ids into a list
        clusters_df['seq_ids'] = clusters_df.seq_ids_raw.apply(lambda r: r.split(','))
        
        otu_ids = list(dedup_info_df.seq.unique())
        rep_otu_ids = list(dedup_info_df.otu_id.unique())
        
        new_dedup_info_df = pd.DataFrame(columns=['otu_id', 'seq'])
        
        for idx, cluster in clusters_df.iterrows():
            #new_dedup_info_df = new_dedup_info_df.append(pd.DataFrame(zip([cluster.otu_id]*len(cluster.seq_ids),cluster.seq_ids), columns=['otu_id','seq']),ignore_index=True)
            cluster_block_df = pd.DataFrame(columns=['otu_id','seq'])
            cluster_block_df.seq=cluster.seq_ids
            cluster_block_df.otu_id = cluster.otu_id
            new_dedup_info_df = pd.concat([new_dedup_info_df, cluster_block_df], ignore_index=True)
        dedup_info_df = new_dedup_info_df
             
    placement_db = sqlite3.connect(args.placements_db)
    
    likelihood_threshold = float(args.likelihood_threshold)
    
    # We can merge dedup_info to seq_data to get the start of our eventual table
    
    seq_data = pd.merge(seq_info_df, dedup_info_df, how='inner', left_on='seq', right_on='seq')
    # Add empty of the remaining columns
    seq_data['name']= ""
    seq_data['ncbi_tax_id']= ""
    seq_data['ncbi_rank']= ""
    
    
    # Now onto the classifications that are located in our database
    
    # First get the rank-order list used in the creation of this set, starting with the shallowest rank at the top 
    rank_order = pd.read_sql("SELECT rank as rank_name, rank_order from ranks order by rank_order DESC", placement_db)
    ranks = [str(r) for r in rank_order.rank_name]
    
    
    for otu_id, otu_block in seq_data.groupby('otu_id'):
        print otu_id
    
        otu_placements = pd.read_sql("""SELECT mcc.want_rank as want_rank, mcc.tax_id as ncbi_tax_id, mcc.rank as ncbi_rank, mcc.likelihood,
                                     taxa.tax_name as name, ranks.rank_order as want_rank_order
                                     FROM multiclass_concat as mcc, taxa, ranks
                                     WHERE mcc.tax_id = taxa.tax_id AND mcc.want_rank = ranks.rank
                                     AND mcc.name = \'%s\'
                                     ORDER BY want_rank_order desc""" % otu_id ,placement_db)
        
        if len(otu_placements) == 0:
            # We didn't get any hits for this otu. Liberalize for matching to the name, in case there has been cruft added to the seqID
            otu_placements = pd.read_sql("""SELECT mcc.want_rank as want_rank, mcc.tax_id as ncbi_tax_id, mcc.rank as ncbi_rank, mcc.likelihood,
                                     taxa.tax_name as name, ranks.rank_order as want_rank_order
                                     FROM multiclass_concat as mcc, taxa, ranks
                                     WHERE mcc.tax_id = taxa.tax_id AND mcc.want_rank = ranks.rank
                                     AND mcc.name LIKE \'%s%%\'
                                     ORDER BY want_rank_order desc""" % otu_id ,placement_db)
            
        # Take our best hit that meets our likelihood threshold
        best_hit = otu_placements[otu_placements.likelihood >= likelihood_threshold].iloc[0]
        
        # For all sequences in this OTU, set the name
        seq_data.loc[seq_data.otu_id==otu_id,'name']=best_hit['name']
        seq_data.loc[seq_data.otu_id==otu_id,'ncbi_tax_id']=best_hit.ncbi_tax_id
        seq_data.loc[seq_data.otu_id==otu_id,'ncbi_rank']=best_hit.ncbi_rank
        
    
    seq_data.to_csv(args.bySeqOutput, index=None, columns=['seq','community','name', 'ncbi_tax_id', 'ncbi_rank', 'otu_id'])
    



if __name__ == '__main__':
    main()