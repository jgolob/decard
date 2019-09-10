#!/usr/bin/env python

import argparse
import csv
import re

CMID_re = re.compile('^SYNCM(?P<cmID>\d+)')


def rename(curName,prefix):
    m = CMID_re.search(curName)
    if m:
        name = prefix+str(m.group('cmID'))
    else:
        name = community
        
    return name

def main():
    args_parser = argparse.ArgumentParser()

    args_parser.add_argument('map', help='CSV file mapping reads')
    args_parser.add_argument('--seq_info','-s', help='Name of seq info output file')
    args_parser.add_argument('--labels','-L', help='Name of labels output file')
    args_parser.add_argument('--rename','-r', help='Relabel to something simpler, with this prefix')
    
    
    
    args = args_parser.parse_args()
    
    with open(args.map) as map_f:
        if args.seq_info:
            seq_info_f = open(args.seq_info,'w')
            seq_info_w = csv.writer(seq_info_f)
        else:
            seq_info_w = None
        
        if args.labels:
            labels_f = open(args.labels,'w')
            labels_w = csv.writer(labels_f)
            labels_w.writerow(['specimen','labels'])
            specimens = set()
        else:
            labels_w = None
        
        if labels_w or seq_info_w:
            reader = csv.DictReader(map_f)
            for row in reader:
                if args.rename:
                    name = rename(row['community'], args.rename)
                else:
                    name = community
                
                if labels_w:
                    if not name in specimens:
                        labels_w.writerow([name,name])
                        specimens.add(name)
                
                if seq_info_w:
                    seq_info_w.writerow([row['seqID'], name])
                    
                
        
        if labels_w:
            labels_f.close()
        
        if seq_info_w:
            seq_info_f.close()
        
        map_f.close()
    
       
        
            
if __name__ == "__main__":
    main()