#!/usr/bin/env python
import argparse
import pandas as pd
import sys
import csv

    
def main():
    args_parser = argparse.ArgumentParser()

    args_parser.add_argument('input', help="OTU table to turn into tallies-wide format")
    args_parser.add_argument('tallieswide', help='Where to put output tallies wide file')
    args_parser.add_argument('--seqinfo', '-si', help='(optional) Headerless CSV file mapping sequence ID to community.')
    args_parser.add_argument('--suffix','-S',help='(optional) Suffix to add to specimen / label to keep track of condition')

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
            
        
    
    tallies_wide = pd.DataFrame(index = otu_table.groupby(['name','ncbi_tax_id','ncbi_rank']).count().index)

    
    for comm_name in otu_table.community.unique():
        tallies_wide[str(comm_name)] = 0
    
    for mi, block in otu_table.groupby(['name','ncbi_tax_id','ncbi_rank']):
        org_counts = block.groupby('community').count().weight
        for comm, count in org_counts.iteritems():
            tallies_wide.loc[mi][comm]=count
    
    if args.suffix:
        suffix = args.suffix
    else:
        suffix=""
    tallies_wide.fillna(0,inplace=True)
    tallies_wide_f = open(args.tallieswide,'w')
    tallies_wide_w = csv.writer(tallies_wide_f)
    # Headers
    
    tallies_wide_w.writerow(['tax_name','tax_id','rank']+[n+suffix for n in tallies_wide.columns])
    tallies_wide_w.writerow(['specimen','','']+[n+suffix for n in tallies_wide.columns])
    tallies_wide_w.writerow(['label','','']+[n+suffix for n in tallies_wide.columns])
    
    for classification, row in tallies_wide.iterrows():
        tallies_wide_w.writerow(list(classification)+list(row))
    
    tallies_wide_f.close()
    

    


if __name__ == '__main__':
    main()