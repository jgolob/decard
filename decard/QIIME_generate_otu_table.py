#!/usr/bin/env python
import argparse
import pandas as pd
import json
import re
import csv
from TaxonCache import TaxonCache
import biom

re_taxonomy_choices = [re.compile('D_(?P<rankID>\d+)__(?P<Name>.+)'),
                       re.compile('(?P<rankID>\w)__(?P<Name>.+)')
                       ]

re_badNames = [
    re.compile('Incertae Sedis', re.IGNORECASE),
    re.compile('uncultured', re.IGNORECASE),
    re.compile('family', re.IGNORECASE),
    re.compile('unclassified', re.IGNORECASE),
]


rankOrders = {
    '0':  'kingdom',
    '1':  'phylum',
    '2':  'class',
    '3':  'order',
    '4':  'family',
    '5':  'genus',
    '6':  'species',
    'k':  'kingdom',
    'p':  'phylum',
    'c':  'class',
    'o':  'order',
    'f':  'family',
    'g':  'genus',
    's':  'species',


}

def get_best_name(row):
    taxonomy_list = list(row.taxonomyList)
    best_name = taxonomy_list.pop()
    best_name = best_name.replace("_"," ")
    # If we're at the species level, some of the taxonomies neglect to include
    # the genus name. Here we'll try to correct that
    if row.ncbi_rank=='species' and len(best_name.split(' ')) < 2 :
        best_name = taxonomy_list[-1]+" "+best_name
        row.taxonomyList[-1]=best_name

    while taxonomy_list and (best_name =='' or len([r for r in re_badNames if r.search(best_name)]) > 0):
        best_name = taxonomy_list.pop()
    
    return best_name


def get_tax_id(start_taxonomy_list,taxons):
    # copy our list, as we'll be popping from it
    tax_id = 1
    taxonomy_list = list(start_taxonomy_list)
    # Strip off less useful things
    best_name = taxonomy_list[-1]
    while taxonomy_list and (best_name =='' or len([r for r in re_badNames if r.search(best_name)]) > 0):
        best_name = taxonomy_list.pop()
    
    # Try to see if our best name has a match
    tax = taxons.lookup_tax_name(best_name)
    if tax:
        return tax['TaxId']
    
    
    while taxonomy_list and tax_id ==1:
        cur_name = taxonomy_list.pop()
        tax = taxons.lookup_tax_name(cur_name)
        if tax:
            tax_id = tax['TaxId']
    
    return tax_id
    


def detect_taxonomy_string_type(rawTaxonomyString):
    
    for re_taxonomy in re_taxonomy_choices:
        if re_taxonomy.search(rawTaxonomyString):
            return re_taxonomy
        
    return None

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

    args_parser.add_argument('input')
    args_parser.add_argument('bySeqOutput')
    args_parser.add_argument('--email','-e',required=True)
    args_parser.add_argument('--seq_info','-s', help="CSV file mapping sequences to specimens")
    args_parser.add_argument('--url','-u',help='Base URL to REST db with taxonomic information. If not provided, we will default to NCBI')

    args = args_parser.parse_args()
    
    if args.url:
        taxons = TaxonCache(source='custom', url = args.url, email=args.email)
    else:
        taxons = TaxonCache(email=args.email)
        
    
    try: # Try first load things in with the janky biom-format.Table class.
        biom_table = biom.load_table(args.input)
        seqs_df = pd.DataFrame(biom_table.ids('sample'), columns=['seq'])
        if len(seqs_df) < 1000:
            print "I'm worried this isn't outputting by sequence. Proceeding without confidence"
            
        otu_df = pd.DataFrame(biom_table.ids('observation'), columns=['otu_id'])
        otu_df['metadata'] =otu_df.otu_id.apply(lambda otu_id: biom_table.metadata(id=otu_id,axis='observation'))
        
        otu_matrix = pd.DataFrame([r for r in biom_table.iter_data(dense=False)], columns=['raw'])
        otu_matrix['otu_index'] = otu_matrix.raw.apply(sparse_matrix_to_otu_index)
        otu_matrix['weight'] = otu_matrix.apply(get_weight, axis=1)
        otu_matrix.drop('raw',1,inplace=True)
        
        matrix_type='sparse'
        
    except: # Next try to load it in as JSON. Allow the error to cascade if needed. 
        in_f = open(args.input)
        biom_raw = json.load(in_f)
        in_f.close()
        seqs_df = pd.DataFrame(biom_raw['columns'])
        seqs_df.rename(columns={'id': 'seq'}, inplace=True)
        # Get rid of the extraneous metadata 
        seqs_df.drop('metadata',1, inplace=True)
        
        otu_df = pd.DataFrame(biom_raw['rows'])
        otu_df.rename(columns={'id': 'otu_id'}, inplace=True)
        otu_matrix = pd.DataFrame(biom_raw['data'])
        matrix_type = biom_raw['matrix_type']
        otu_matrix.columns = ['otu_index','seq','weight']
        

    

    # fix up our taxonomy strings and get a name and NCBI
    # break up and clean up the taxonomy strings
    # First figure out which kind of taxonomy string we have by searching a bit
    re_taxonomy = detect_taxonomy_string_type(otu_df.iloc[0].metadata['taxonomy'][0])
    
    otu_df['taxonomyList'] = otu_df.metadata.apply(lambda metadata: [re_taxonomy.search(raw_tax).group('Name') for raw_tax in metadata['taxonomy'] if re_taxonomy.search(raw_tax)] )
    otu_df['ncbi_rank'] = otu_df.metadata.apply(lambda metadata: rankOrders.get([re_taxonomy.search(raw_tax).group('rankID') for raw_tax in metadata['taxonomy'] if re_taxonomy.search(raw_tax)][-1]))
    otu_df['name'] = otu_df.apply(get_best_name, axis=1)
    otu_df['ncbi_tax_id'] = otu_df.taxonomyList.apply(lambda tid: get_tax_id(tid, taxons))
    otu_df['taxonomy_string'] = otu_df.taxonomyList.apply(lambda tl: '; '.join(tl))
    otu_df.drop('metadata',1,inplace=True)
    
    
    
    
    if matrix_type =='sparse':
        
        # Sequences are indexed by row, so we can do this the easy way
        otu_matrix['seq'] = seqs_df.seq
        # for otus, we need to do a merge, using the index on the otu dataframe to the value in the OTU id column
    
        otu_matrix = pd.merge(otu_matrix, otu_df, how='inner', left_on='otu_index', right_index=True)
        
    if args.seq_info:
        seq_info_df = pd.read_csv(args.seq_info, header=0, names=['seq','community'])
        otu_matrix= pd.merge(otu_matrix,seq_info_df, how='inner', left_on='seq', right_on='seq')
    else:
        otu_matrix['community'] = ""
        
    otu_matrix.drop('otu_index', inplace=True, axis=1)
    otu_matrix.drop('taxonomyList', inplace=True, axis=1)
    otu_matrix.to_csv(args.bySeqOutput, index=None)


if __name__ == '__main__':
    main()