#!/usr/bin/env python

import argparse
import csv
import re

CMID_re = re.compile('^SYNCM(?P<cmID>\d+)')

def rename(curName, prefix):
    m = CMID_re.search(curName)
    if m:
        name = prefix+str(m.group('cmID'))
    else:
        name = community
        
    return name

def main():
    args_parser = argparse.ArgumentParser()

    args_parser.add_argument('map', help='CSV file mapping reads')
    args_parser.add_argument('output', help='Where to put output CSV')
    args_parser.add_argument('--rename','-r', help='Relabel to something simpler, with this prefix')
    args_parser.add_argument('--suffix','-s', help='Suffix to add to name')
    
    args = args_parser.parse_args()
    
    if args.suffix:
        suffix = args.suffix
    else:
        suffix = ''
    
    map_data = []
    with open(args.map) as in_f:
        reader = csv.DictReader(in_f)
        for row in reader:
            map_data.append(row)
        
        in_f.close()
    
    organisms = {}
    communities = {}
    
    for data in map_data:
        community = data['community']
        
        if args.rename:
            community = rename(community, args.rename)
            
        community = community+suffix
        
        taxID = data['ncbi_tax_id']
        organisms[taxID]=data['organism']
        if not community in communities:
            communities[community]={}
        
        try:
            communities[community][taxID]+=1
        except KeyError:
            communities[community][taxID]=1
        
    
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