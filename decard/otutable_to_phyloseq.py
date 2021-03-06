#!/usr/bin/env python
"""
        Goal: Take an otu table as input, and covert it to something that can be dragged through the R
        Phyloseq package.
        Input: xxxx.otutable.csv
        
"""
import argparse
import pandas as pd

def get_species(row):
    """
        Try to figure out if we have a species name
    """
    if len(row.taxonomy_list)> 6: # If our taxonomy is long enough, use the seventh position as our species name
        return row.taxonomy_list[6]
    elif row.organism != row.taxonomy_list[-1]: # Do we have an organism name different than our last spot in the taxonomy list? If so, use that
        return row.organism
    else:
        return "" # Nothing useful. 
    

def main():
        args_parser = argparse.ArgumentParser()
        
        
        args_parser.add_argument('input', help="OTU table in csv format")
        args_parser.add_argument('output_prefix', help="OUTPUT prefix to use (including path) for the two tables")
        args_parser.add_argument('--seqinfo', '-si', help='(optional) Headerless CSV file mapping sequence ID to community.')
        
        args = args_parser.parse_args()

            
        # Open our OTU table
    
        otu_table = pd.read_csv(args.input)
        # Do a bit of normalizing here, first make sure that our OTU table has a weight column. If not, default to a weight of 1
        if not 'weight' in otu_table.columns:
            otu_table['weight'] = 1
            
        # Next be sure that we have a community column, and that it actually has community data in it.
        if not 'community' in otu_table.columns or len(otu_table.community[otu_table.community.isnull()]) > 0:
            # if we do not have seq info data, stop here.
            if not args.seqinfo:
                sys.stderr.write("No community information avaiable. Consider providing a seqinfo file")
                return -1
            else:
                # Load the seq info
                seq_info_df = pd.read_csv(args.seqinfo, names=['seq','community'])
                # Drop the bad community column
                try:
                    otu_table.drop('community', axis=1, inplace=True)
                except ValueError: # Wasn't one in the first place. Fine
                    pass
                otu_table = pd.merge(otu_table,seq_info_df,on='seq')
            
        #taxons = TaxonCache(source='custom', email='jgolob@fhcrc.org', url='http://fredricks-rome.fhcrc.org/taxon/JSON/ncbi/')
        
        # Create the OTU matrix, with each community getting a column, each OTU getting a row.
        otu_matrix_df = pd.DataFrame(columns=otu_table.community.unique(), index=otu_table.otu_id.unique())
        
        # Fill our matrix by looping through, first grouping by OTU_ID
        for otu_id, otu_block in otu_table.groupby('otu_id'):
            # Then get the total weight for each OTU, by community
            community_weights = otu_block.groupby('community').sum().weight
            for comm, weight in community_weights.iteritems():
                # Use the loc accessor to put in the value
                otu_matrix_df.loc[otu_id,comm]= weight
        
        # Fill the remaining with zeros
        otu_matrix_df.fillna(0, inplace=True)
        
        # Now work on our tax table, where the index is the otu_id.
        
        # first get our uniqe otu_id / taxon strings
        otu_strings = otu_table.drop_duplicates(['otu_id','taxonomy_string','name'])
        
        # Next set up our taxon_df
        taxon_df = pd.DataFrame()
        taxon_df['otu_id'] = otu_strings.otu_id
        taxon_df['taxonomy_string'] = otu_strings.taxonomy_string
        taxon_df['organism'] = otu_strings['name']
        # Then make our OTU ID the index
        taxon_df.set_index('otu_id', inplace=True)
         
        # Split our strings into lists
        taxon_df['taxonomy_list'] = taxon_df.taxonomy_string.apply(lambda ts: ts.split('; '))
        
        taxon_df['Domain'] = taxon_df.taxonomy_list.apply(lambda tl: tl[0] if len(tl) > 0 else "")
        taxon_df['Phylum'] = taxon_df.taxonomy_list.apply(lambda tl: tl[1] if len(tl) > 1 else "")
        taxon_df['Class'] = taxon_df.taxonomy_list.apply(lambda tl: tl[2] if len(tl) > 2 else "")
        taxon_df['Order'] = taxon_df.taxonomy_list.apply(lambda tl: tl[3] if len(tl) > 3 else "")
        taxon_df['Family'] = taxon_df.taxonomy_list.apply(lambda tl: tl[4] if len(tl) > 4 else "")
        taxon_df['Genus'] = taxon_df.taxonomy_list.apply(lambda tl: tl[5] if len(tl) > 5 else "")
        taxon_df['Species'] = taxon_df.apply(get_species,axis=1)
        
                
        
        # OUTPUT TIME!
        taxon_df.to_csv(args.output_prefix+'.tax_table.csv', columns=['Domain', 'Phylum', 'Class', 'Order', 'Family','Genus', 'Species'])
        otu_matrix_df.to_csv(args.output_prefix+".otu_matrix.csv")

if __name__ == '__main__':
    main()
