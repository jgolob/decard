#!/usr/bin/env python
"""
    This module takes a map file, and checks to be sure the ncbi_tax_id for a given row is actually correct. It's a QC step.
    INPUT: A map file
    
"""


import argparse
import csv
from TaxonCache import TaxonCache


def main():
    args_parser = argparse.ArgumentParser()

    args_parser.add_argument('map_file',  help='Sequences map file that identifies each seq')
    args_parser.add_argument('--email','-e',required=True)
    args_parser.add_argument('--url','-u',help='Base URL to REST db with taxonomic information. If not provided, we will default to NCBI')
    
    
    args = args_parser.parse_args()
    
    if args.url:
        taxons = TaxonCache(source='custom', url = args.url, email=args.email)
    else:
        taxons = TaxonCache(email=args.email)
    
    with open(args.map_file) as map_f:
        map_r = csv.DictReader(map_f)
        with open(args.labels,'w') as labels_f:
            taxons
            
            
            labels_f.close()
        map_f.close()
            
if __name__ == "__main__":
    main()
