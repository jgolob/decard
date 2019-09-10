#!/usr/bin/env python
"""
        Goal: Take an otu table as input, and covert it to something that can be dragged through the R
        Phyloseq package.
        Input: xxxx.otutable.csv
        
"""
import argparse
import pandas as pd
import csv
from TaxonCache import TaxonCache
import collections


def flatten(l):
    for el in l:
        if isinstance(el, collections.Iterable) and not isinstance(el, basestring):
            for sub in flatten(el):
                yield sub
        else:
            yield el

# ---
def lookup_row_at_rank(taxons, row, rank):
        ncbi_ids = row.otu_id.split(',')
        lineages = [taxons.lookup_tax_id(tid)['LineageEx'] for tid in ncbi_ids]
        return ';'.join(set(flatten([[tax['ScientificName'] for tax in lineage if tax["Rank"]==rank] for lineage in lineages])))
        


def main():
        args_parser = argparse.ArgumentParser()
        
        args_parser.add_argument('tallies_wide', help='INPUT  tallies_wide file, in CSV format')
        args_parser.add_argument('output_prefix', help="OUTPUT prefix to use (including path) for the two tables")
        args_parser.add_argument('--email','-e',required=True)
        args_parser.add_argument('--url','-u',help='Base URL to REST db with taxonomic information. If not provided, we will default to NCBI')

        args = args_parser.parse_args()
        # Load the taxon cache
        
        if args.url:
            taxons = TaxonCache(source='custom', url = args.url, email=args.email)
        else:
            taxons = TaxonCache(email=args.email)
        
        
        #taxons = TaxonCache(source='custom', email='jgolob@fhcrc.org', url='http://fredricks-rome.fhcrc.org/taxon/JSON/ncbi/')
        
        
        # Load data
        tw_f = open(args.tallies_wide)
        tw_reader = csv.reader(tw_f)
        header_1 = tw_reader.next()
        header_2 = tw_reader.next()
        header_3 = tw_reader.next()
        labels = header_3[3:]
        
        taxon_table = []
        otu_table = [['']+labels]
        for row in tw_reader:
                otu_id = row[1]
                counts = row[3:]
                taxon_table.append(row[:3])
                otu_table.append([otu_id]+counts)
        
        otu_df = pd.DataFrame(otu_table[1:], columns=['otu_id']+otu_table[0][1:])
        otu_df.set_index('otu_id', inplace=True)
        
        # Now work on our tax table
        taxon_df = pd.DataFrame(taxon_table, columns=['organism', 'otu_id', 'rank'])
        taxon_df['Domain'] = taxon_df.apply(lambda row: lookup_row_at_rank(taxons,row,'superkingdom'), axis=1)
        taxon_df['Phylum'] = taxon_df.apply(lambda row: lookup_row_at_rank(taxons,row,'phylum'), axis=1)
        taxon_df['Class'] = taxon_df.apply(lambda row: lookup_row_at_rank(taxons,row,'class'), axis=1)
        taxon_df['Order'] = taxon_df.apply(lambda row: lookup_row_at_rank(taxons,row,'order'), axis=1)
        taxon_df['Family'] = taxon_df.apply(lambda row: lookup_row_at_rank(taxons,row,'family'), axis=1)
        taxon_df['Genus'] = taxon_df.apply(lambda row: lookup_row_at_rank(taxons,row,'genus'), axis=1)
        taxon_df['Species'] = taxon_df.organism
        taxon_df.set_index('otu_id', inplace=True)
                
        
        # OUTPUT TIME!
        taxon_df.to_csv(args.output_prefix+'.tax_table.csv', columns=['Domain', 'Phylum', 'Class', 'Order', 'Family','Genus', 'Species'])
        otu_df.to_csv(args.output_prefix+".otu_matrix.csv")

if __name__ == '__main__':
    main()
