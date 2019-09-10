#!/usr/bin/env python

import argparse

import gzip
import re
from Bio import SeqIO
import os
import sys
from TaxonCache import TaxonCache


re_silva_description = re.compile('(?P<accession>[\w]+)\.(?P<acc_loc>[\d\.]+) (?P<taxonomy>[\w; ]+)')

unambiguousBases = ['a','c','t','g','u']


def no_ambiguous_bases(seq):
    return (  len([b for b in seq.lower() if not b in unambiguousBases]) == 0)

def main():
    args_parser = argparse.ArgumentParser()

    args_parser.add_argument('--directory', '-d', help='Directory where to put genus directory', required=True)
    args_parser.add_argument('--genus','-g', help='Genus from which to select sequences', required=True)
    args_parser.add_argument('--reference','-r', help='SILVA reference file', required=True)
    args_parser.add_argument('--mock', '-m', help='Mock run. Do not modify the FASTA file and limit how many records we go after', action='store_true')
    args_parser.add_argument('--email','-e',required=True)
    args_parser.add_argument('--url','-u',help='Base URL to REST db with taxonomic information. If not provided, we will default to NCBI')
    
    args = args_parser.parse_args()

    silva_f = gzip.open(args.reference)
    silva_sr = SeqIO.parse(silva_f,'fasta')
    
    target_genus = args.genus
    target_genus = target_genus.lower()
    
    if args.url:
        taxons = TaxonCache(source='custom', url = args.url, email=args.email)
    else:
        taxons = TaxonCache(email=args.email)
    
    out_f_handles = {}
    
    for sr in silva_sr:
        m = re_silva_description.search(sr.description)
        if m:
            taxonomyString = m.group('taxonomy')
            taxonomyList = taxonomyString.split(';')
            genus = taxonomyList[-2]
            if genus.lower() == target_genus.lower():
                if not taxons.is_bad_name(taxonomyList[-1]) and no_ambiguous_bases(sr.seq):
                    organismName = taxons.get_best_name(taxonomyList)
                    organism = " ".join(organismName.split(' ')[:2]) # Get rid of strain junk for our filename
                    print organism
                    tax = taxons.lookup_taxonomy(taxonomyList)
                    taxID = tax['TaxId']
                    # Make description conform to expected
                    sr.description='description=\"'+sr.description+'\" organism=\"'+organism+'\" taxonomy=\"'+';'.join(taxonomyList[:-1])+'\" ncbi_tax_id=\"'+str(taxID)+'\"'
                    
                    if not args.mock:
                        try:
                            of = out_f_handles[organism]
                        except KeyError:
                            if not os.path.isdir(args.directory+'/'+args.genus):
                                os.mkdir(args.directory+'/'+args.genus)
                            of = open(args.directory+'/'+args.genus+'/'+organism+".fasta",'a')
                            out_f_handles[organism] = of
                        
                        SeqIO.write(sr,of,'fasta')
                        
                    else:
                        sys.stdout.write(str(sr.description)+'\n')
    
    
    for org in out_f_handles:
        out_f_handles[org].close()



if __name__ == "__main__":
    main()