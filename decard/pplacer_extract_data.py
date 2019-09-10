#!/usr/bin/env python
"""
    An alternative method to extract classifications from the placement db, using a somewhat fancier method.
    
    The relevant bits of the schema:
    -- Multiclass is where the classifications live.
        -- name is the cluster's representitive sequence ID
        -- want_rank specifies the goal rank for the classifer
        -- rank is the actual rank of the classified taxon
        -- tax_id points to the taxa table, giving us the name. By default this is an NCBI taxon id and externally valid as well (int)
        -- likelihood is the likelihood this classification being true
        -- placement_id links to the placements table, and gives you by which method the classification was completed, and which run
        -- A given cluster sequence thus will have multiple classifications for different want_ranks, and potentially for classifcation methods and runs
        
    CREATE TABLE multiclass (
            placement_id INTEGER REFERENCES placements (placement_id) NOT NULL,
            name TEXT NOT NULL,
            want_rank TEXT REFERENCES ranks (rank) NOT NULL,
            rank TEXT REFERENCES ranks (rank) NOT NULL,
            tax_id TEXT REFERENCES taxa (tax_id),
            likelihood REAL NOT NULL
          );
        

    CREATE TABLE multiclass_concat(
            placement_id INT,
            name TEXT,
            want_rank TEXT,
            tax_id,
            rank,
            likelihood,
            id_count
          );
    
    -- Placement names is how one can get the mass (how many sequences this centroid sequence is ultimately representing)
    CREATE TABLE placement_names (
                   placement_id INTEGER NOT NULL,
                   name TEXT NOT NULL,
                   mass REAL NOT NULL,
                   PRIMARY KEY (name));
                   
    -- As described in the multiclass table. Specifies which classifier and run
    
    CREATE TABLE placements (
            placement_id INTEGER PRIMARY KEY AUTOINCREMENT,
            classifier TEXT NOT NULL,
            run_id INTEGER REFERENCES runs (run_id) NOT NULL
          );
    
    -- Ranks in our taxonomy, along with their rank-order
    CREATE TABLE ranks (
            rank TEXT PRIMARY KEY NOT NULL,
            rank_order INTEGER
          );
          
    -- A flow in of the original sequence IDs and from which sequence each came.
        -- Sadly, we do not have a real deduplication / derep table that specifies to which cluster these sequences were placed. 
    CREATE TABLE seq_info (
        name VARCHAR(37) NOT NULL, 
        specimen VARCHAR(10) NOT NULL
    );
    
    -- A subset of the NCBI taxonomy. 
    CREATE TABLE taxa (
            tax_id TEXT PRIMARY KEY NOT NULL,
            tax_name TEXT NOT NULL,
            rank TEXT REFERENCES ranks (rank) NOT NULL
          );
    
    
    Our method:
        
        For each centroid sequence:
            Start with the lowest rank (deepest into the taxonomic tree) for want_rank, and step upwards until
                the likelihood of classification at this rank is above a specified threshold
                If multiple hits at this rank, see if they are similar
                    -> If not, step up
                    -> If yes, consider taking the parent ranking
                    
            
            Have the cetroid classified to this sequence. 
        
"""


import argparse
import sqlite3
import pandas as pd
import csv
import numpy as np

def main():
    args_parser = argparse.ArgumentParser()

    args_parser.add_argument('placementDB', help="Placement DB")
    args_parser.add_argument('--threshold','-t', default=0.99, help='Threshold likelihood to call')
    args_parser.add_argument('--by_sequence_output','-S', help='CSV filename into which to put the by-species classification')
    args_parser.add_argument('--dedup_info','-D', help='Deduplication info')
    
    args = args_parser.parse_args()
    
    threshold = float(args.threshold)
       
    if True:
        with sqlite3.connect(args.placementDB) as conn:
            
            # Get the ranks
            print "Loading Ranks"
            rank_df = pd.read_sql('SELECT rank, rank_order from ranks order by rank_order asc',conn, index_col = 'rank')
            ranks = list(rank_df.index)
            print "\t Ranks loaded"
        
            
            # Get the centroid IDs
            
            print "Loading centroid data"
            centroid_ids_df = pd.read_sql('SELECT DISTINCT name from placement_names', conn)
            print "\t centroid data loaded"
            
            print "Loading placements (this can take a while)"
            placements=pd.read_sql(""" 
                    SELECT mcc.name as clusterID, mcc.rank as rank, mcc.likelihood as likelihood, mcc.want_rank as want_rank,
                        taxa.tax_name as ncbi_tax_name, mcc.tax_id as ncbi_tax_id
                        FROM multiclass_concat as mcc JOIN taxa ON mcc.tax_id = taxa.tax_id
                        WHERE mcc.likelihood >= %f""" % (threshold),conn)
            print '\t Done loading placements (whew)'

            
            # Initialize our variables to mak ethis work
            classified = pd.DataFrame()
            classified['clusterID']=[]
            not_classified = centroid_ids_df.name
                        
            while ranks and len(not_classified)>0:
                rank = ranks.pop()
                print len(classified), 'of', len(centroid_ids_df), ' reads classified. Attempting ', rank
                
                placements_at_rank = placements.query('want_rank == @rank')
                
                
                duplicated = placements_at_rank['clusterID'].duplicated()
                duplicated_at_rank = placements_at_rank[duplicated]
                if len(duplicated_at_rank) >0:
                    print "DUPLICATES!"
                # deal with duplicated here...
                unique_at_rank = placements_at_rank[~duplicated]
                newly_captured_ids = np.setdiff1d(unique_at_rank.clusterID, classified.clusterID)
                
                # Now selectively add the newly captured IDs
                classified = pd.concat([classified,unique_at_rank.query('clusterID in @newly_captured_ids')])        
                
                
                not_classified = np.setdiff1d(centroid_ids_df.name,unique_at_rank.clusterID)
                print '\t', len(newly_captured_ids), ' newly classified. ', len(not_classified), ' to go'
            
            assert(len(classified)==len(centroid_ids_df))
            
            classified['organism_name']=classified['ncbi_tax_name']
            
            
            """
                For classifications of a rank above species, ideally we would break these into individual names (e.g. novel clostridial species #2)
                
                There is a bit of quirkiness to how the pplacer outputs OTUs. In this case, it's an OTU PER specimen that we've actually classified. 
            """
            
            
            if args.dedup_info:
                dedup_info_df = pd.read_csv(args.dedup_info, header=None, names=['centroid', 'source', 'mass'])
                classified_dedup = classified.merge(dedup_info_df, left_on='clusterID', right_on='source')
                
                rank_order = rank_df.rank_order
                species_rank_order = int(rank_order['species'])
                
                not_species = classified_dedup[list(rank_order[classified['rank']] < species_rank_order)]
        
                for (tax_id, tax_name, rank), block in not_species.groupby(['ncbi_tax_id','ncbi_tax_name', 'rank']):
                    print tax_name, len(block.centroid.unique())
                    # for each member of this block, generate an organism name and assign it. 
                    for i, row in enumerate(block.iterrows()):
                        classified.iloc[row[0]].organism_name = tax_name+" us"+unicode(i)

            
            
            
            # Great, we've now gotten our best possible classifications. Extract the unique ncbi_tax_id and ncbi_taxon_name rank triplets
            taxons = classified.groupby(['ncbi_tax_id','ncbi_tax_name', 'rank']).count()
            
                
            

if __name__ == '__main__':
    main()