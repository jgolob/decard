#!/usr/bin/env python
"""
        Goal: Take an OTU table as input, and add the taxon_string column
        Input: xxxx.otutable.csv
        Output: yyyy.otutable.csv
        
"""
import argparse
import pandas as pd
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
        try:
            ncbi_ids = row.ncbi_tax_id.split(',')
        except AttributeError:
            ncbi_ids = [str(row.ncbi_tax_id)]
        lineages = [taxons.lookup_tax_id(tid)['LineageEx'] for tid in ncbi_ids]
        return ' / '.join(set(flatten([[tax['ScientificName'] for tax in lineage if tax["Rank"]==rank] for lineage in lineages])))
        


def main():
        args_parser = argparse.ArgumentParser()
        
        args_parser.add_argument('input', help='INPUT  otutable file, in CSV format')
        args_parser.add_argument('output', help="OUTPUT modified otutable in CSV format, with added tax string")
        args_parser.add_argument('--email','-e',required=True)
        args_parser.add_argument('--url','-u',help='Base URL to REST db with taxonomic information. If not provided, we will default to NCBI')

        args = args_parser.parse_args()
        # Load the taxon cache
        
        if args.url:
            taxons = TaxonCache(source='custom', url = args.url, email=args.email)
        else:
            taxons = TaxonCache(email=args.email)
        
        
        #taxons = TaxonCache(source='custom', email='jgolob@fhcrc.org', url='http://fredricks-rome.fhcrc.org/taxon/JSON/ncbi/')
        
        
        # Open our OTU table
    
        otu_table = pd.read_csv(args.input)
        
        
        otu_table['Domain'] = otu_table.apply(lambda row: lookup_row_at_rank(taxons,row,'superkingdom'), axis=1)
        otu_table['Phylum'] = otu_table.apply(lambda row: lookup_row_at_rank(taxons,row,'phylum'), axis=1)
        otu_table['Class'] = otu_table.apply(lambda row: lookup_row_at_rank(taxons,row,'class'), axis=1)
        otu_table['Order'] = otu_table.apply(lambda row: lookup_row_at_rank(taxons,row,'order'), axis=1)
        otu_table['Family'] = otu_table.apply(lambda row: lookup_row_at_rank(taxons,row,'family'), axis=1)
        otu_table['Genus'] = otu_table.apply(lambda row: lookup_row_at_rank(taxons,row,'genus'), axis=1)
        
        otu_table['taxonomy_string'] = otu_table.apply(lambda r: '; '.join([t for t in [
                r.Domain,
                r.Phylum,
                r.Class,
                r.Order,
                r.Family,
                r.Genus
            ] if t]), axis=1)                
        
        # Catch anything that did't yield anything
        otu_table.taxonomy_string.loc[otu_table.taxonomy_string.isnull()]=otu_table[otu_table.taxonomy_string.isnull()]['name']
        
        
        # Cleanup
        otu_table.drop('Domain', inplace=True, axis=1)
        otu_table.drop('Phylum', inplace=True, axis=1)
        otu_table.drop('Class', inplace=True, axis=1)
        otu_table.drop('Order', inplace=True, axis=1)
        otu_table.drop('Family', inplace=True, axis=1)
        otu_table.drop('Genus', inplace=True, axis=1)
        
        # OUTPUT TIME!
        otu_table.to_csv(args.output, index=None)

if __name__ == '__main__':
    main()
