#!/usr/bin/env python
import argparse
import pandas as pd
from Bio import SeqIO
import json
import re
import csv
from TaxonCache import TaxonCache

re_taxonomy_choices = [re.compile('D_(?P<rankID>\d+)__(?P<Name>.+)'),
                       re.compile('(?P<rankID>\w)__(?P<Name>.+)')
                       ]

re_certainty = re.compile('(?P<name>.+)\((?P<certainty>\d+)\)')

rankOrders = {
    '1':  'kingdom',
    '2':  'phylum',
    '3':  'class',
    '4':  'order',
    '5':  'family',
    '6':  'genus',
    '7':  'species',
    'k':  'kingdom',
    'p':  'phylum',
    'c':  'class',
    'o':  'order',
    'f':  'family',
    'g':  'genus',
    's':  'species',


}



def detect_taxonomy_string_type(rawTaxonomyString):
    
    for re_taxonomy in re_taxonomy_choices:
        if re_taxonomy.search(rawTaxonomyString):
            return re_taxonomy
        
    return None


def taxID_from_tax(tax):
    if tax:
        return int(tax['TaxId'])
    else:
        return 1

def load_classifications(classification_fn, taxons):
    """
        Input: the filename of our final_taxonomy.
            This is a tab delimited CSV-like file
            Col 1 is the OTU ID
            Col 2 count (weight)
            Col 3 is the proposed taxonomy for this OTU, in the format:
            KINGDOM(likihood); PHYLUM(likelihood);.....
            eg.
            Bacteria(100);Firmicutes(100);Clostridia(100);Clostridiales(100);Ruminococcaceae(100);Gemmiger(93);
            
            We will read in this file, and then churn through this raw string to break it down into something closer to our standard columns:
            TaxonomyList -> A list, with each rank in a seperate block
            ncbi_rank -> Rank of classification (Genus in our example)
            taxonomy_string -> Semi-delimited string ala
            Bacteria ;Firmicutes ;Clostridia ;Clostridiales ;Ruminococcaceae ;Gemmiger
            name -> A best name, filtering out nonsense
                "Gemminger"
            ncbi_tax_id
    """
    classification_df = pd.read_csv(classification_fn, header=0, sep='\t',names=['otu_id','weight','taxonomy_raw'])
    classification_df['taxonomyListRaw'] = classification_df.taxonomy_raw.apply(lambda ts: ts.split(';'))
    classification_df['taxonomyList'] = classification_df.taxonomyListRaw.apply(lambda tl: [re_certainty.search(t).group('name') for t in tl if re_certainty.search(t)])
    classification_df['ncbi_rank'] = classification_df.taxonomyList.apply(lambda tl: rankOrders.get(str(len(tl))))
    classification_df['taxonomy_string'] = classification_df.taxonomyList.apply(lambda tl: '; '.join(tl))
    classification_df['name'] = classification_df.taxonomyList.apply(lambda tl: taxons.get_best_name(tl))
    classification_df['ncbi_tax_id'] = classification_df.taxonomyList.apply(lambda tl: taxID_from_tax(taxons.lookup_taxonomy(tl)) )
    classification_df.set_index('otu_id', inplace=True)
    
    
    return classification_df

def load_otu_seqs(list_fn, target_percent=0.03):
    list_f = open(list_fn)
    header = list_f.readline()
    header = header.replace('\n','')
    header = header.split('\t')
    otu_names = header[2:]
    for line in list_f:
        line = line.replace('\n','')
        line = line.split('\t')
        target = line[0]
        if target == str(target_percent):
            num_otu = int(line[1])
            otus = line[2:]
            otu_seqs = [o.split(',') for o in otus]
            list_f.close()
            return otu_seqs

    # If we're here, we weren't able to get at our target percent. Curious. Warn the user and then take ANY data.
    print "Could not find any data at the target percent ", target_percent, "Returning what I can find"
    list_f.seek(0)
    header = list_f.readline()
    try:
        line = list_f.next()
        line = line.replace('\n','')
        line = line.split('\t')
        target = line[0]
        num_otu = int(line[1])
        otus = line[2:]
        otu_seqs = [o.split(',') for o in otus]
    except:
        otu_seqs = []
    list_f.close()
    return otu_seqs

def seqs_for_representitive(otu_id, otu_seqs):
    """
        Like many classification systems, one ends up with a classification for a representitive sequence for an OTU.
        Given seq_ids grouped by list into OTUs, and an OTU ID that is one of the sequences in an OTU, return a list
        of sequences that should be with this classification
    """
    
    for OTU_seq in otu_seqs:
        if otu_id in OTU_seq:
            return OTU_seq
    
    return []

def classification_for_seqs(seqs, classification_df):
    """
        Like many classification systems, one ends up with a classification for a representitive sequence for an OTU.
        Given a list of sequences in an OTU (seqs), pick a classification from a classification dataframe
        Assumes index of classification_df is the seq_id
    """
    for seq_id in seqs:
        try:
            class_idx = classification_df.index.get_loc(seq_id)
            return classification_df.iloc[class_idx]
        except KeyError: # Not the representive
            pass
    
    # Couldn't find anything
    return None
      
def main():
    args_parser = argparse.ArgumentParser()

    args_parser.add_argument('otu_taxonomy', help='OTU->Taxonomy .good.unique.filter.unique.precluster.pick.pick.an.x.cons.taxonomy')
    args_parser.add_argument('dedup_fasta', help='Fasta file that contains all the deduplication info plus sequences (.good.unique.filter.unique.precluster.pick.pick.an.x.fasta)')
    #args_parser.add_argument('final_taxonomy', help='File of representitive sequence and its classification')
    args_parser.add_argument('otu_table',help='output file into which to put the otu table')
    args_parser.add_argument('--email','-e',required=True)
    args_parser.add_argument('--seq_info','-s', help="CSV file mapping sequences to specimens")
    args_parser.add_argument('--url','-u',help='Base URL to REST db with taxonomic information. If not provided, we will default to NCBI')

    args = args_parser.parse_args()
    if args.url:
        taxons = TaxonCache(source='custom', url = args.url, email=args.email)
    else:
        taxons = TaxonCache(email=args.email)
    

    # Load in the dedup info from the fasta file
    dedup_seqs = SeqIO.parse(args.dedup_fasta,'fasta')
    # First column will be seq_ids, second is OTU represented
    dedup_info = [sr.description.split() for sr in dedup_seqs]
    dedup_df = pd.DataFrame(dedup_info, columns=['seq','otu_id'])
    
    
    # Then load in the classifications
    classification_df = load_classifications(args.otu_taxonomy, taxons)
    
    # merge the two to make our nascent OTU table
    otu_df = pd.merge(dedup_df, classification_df, left_on='otu_id', right_index=True, how='left')
    
    seq_info = {}
    if args.seq_info:
        # If we have a mapping of sequence to community, load it here
        seq_info_df = pd.read_csv(args.seq_info, names=['seq','community'], header=-1)
        
        # and merge it
        otu_df = pd.merge(otu_df, seq_info_df, on='seq', how='left')
    else: # Add the community column but leave it empty
        otu_df['community'] = ''

    # And then write out to the csv file
    otu_df.to_csv(args.otu_table, index=False, columns=['seq','community','otu_id', 'name', 'ncbi_rank','ncbi_tax_id', 'taxonomy_string'])
    
    



if __name__ == '__main__':
    main()