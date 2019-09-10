#!/usr/bin/env python

import argparse
import csv
import re
from Bio import SeqIO

CMID_re = re.compile('^SYNCM(?P<cmID>\d+)')
header_re = re.compile('description=\"(?P<description>[^\"]+)\" organism=\"(?P<organism>[^\"]+)\" taxonomy=\"(?P<taxonomy>[^\"]+)\" ncbi_tax_id=\"(?P<ncbi_tax_id>[^\"]+)\"')


def rename(curName, prefix):
    m = CMID_re.search(curName)
    if m:
        name = prefix+str(m.group('cmID'))
    else:
        name = community
        
    return name

def load_targets(target_fn):
    target_data = []
    with open(target_fn) as in_f:
        reader = csv.DictReader(in_f)
        for row in reader:
            target_data.append(row)
        
        in_f.close()
    return target_data

def add_counts(target_data):
    weights = [float(t['weight']) for t in target_data]
    min_weight = min(weights)
    for target in target_data:
        target['count'] = int(round(float(target['weight'])/min_weight))
        

def add_tax_id(target_data):
    for target in target_data:
        seq_recs = SeqIO.parse(target['source_file'],'fasta')
        for sr in seq_recs:
            if sr.id == target['sequence_id']:
                m = header_re.search(sr.description)
                if m:
                    target['ncbi_tax_id'] = m.group('ncbi_tax_id')
                    target['species'] = m.group('organism')
                else:
                    target['ncbi_tax_id'] = 1
                break
        if not 'ncbi_tax_id' in target:
            target['ncbi_tax_id'] = 1
                    
                
            


        

def main():
    args_parser = argparse.ArgumentParser()

    args_parser.add_argument('target', help='CSV file targets for communities')
    args_parser.add_argument('output', help='Where to put output CSV')
    args_parser.add_argument('--rename','-r', help='Relabel to something simpler, with this prefix')
    args_parser.add_argument('--suffix','-s', help='Suffix to add to name')
    
    args = args_parser.parse_args()
    
    if args.suffix:
        suffix = args.suffix
    else:
        suffix = ''
    
    target_data = load_targets(args.target)
    
    # Does this by reference, modifying the original...
    add_counts(target_data)
    
    # Does this by reference, modifying the original...
    add_tax_id(target_data)
    
    
    organisms = {}
    communities = {}
    
    for data in target_data:
        community = data['community_name']
        
        if args.rename:
            community = rename(community, args.rename)
            
        community=community+suffix
        
        taxID = data['ncbi_tax_id']
        organisms[taxID]=data['species']
        if not community in communities:
            communities[community]={}
        
        try:
            communities[community][taxID]+=data['count']
        except KeyError:
            communities[community][taxID]=data['count']
        
    
    with open(args.output,'w') as out_f:
        writer = csv.writer(out_f)
        # Header line 1
        writer.writerow(["tax_name",'tax_id','rank']+[c for c in communities])
        # Header lines 2 and 3
        writer.writerow(["specimen",'','']+[c for c in communities])
        writer.writerow(["label",'','']+[c for c in communities])
        for taxID in organisms:
            row = [organisms[taxID], taxID, 'species']
            for comm in communities:
                row.append(communities[comm].get(taxID,0))
            writer.writerow(row)

        out_f.close()            
    
       
        
            
if __name__ == "__main__":
    main()