#!/usr/bin/env python

import argparse
from Bio import SeqIO
import re
import os
import sys

header_re = re.compile('description=\"(?P<description>[^\"]+)\" organism=\"(?P<organism>[^\"]+)\" taxonomy=\"(?P<taxonomy>[^\"]+)\" ncbi_tax_id=\"(?P<ncbi_tax_id>[^\"]+)\"')

def main():
    args_parser = argparse.ArgumentParser()

    args_parser.add_argument('--directory', '-d', help='Directory where to put genus directory', required=True)
    args_parser.add_argument('--genus','-g', help='Genus from which to select sequences', required=True)
    args_parser.add_argument('--reference','-r', help='FASTA reference file', required=True)
    args_parser.add_argument('--mock', '-m', help='Mock run. Do not modify the FASTA file and limit how many records we go after', action='store_true')
    
    args = args_parser.parse_args()
    
    """
    if not os.path.isdir(args.directory):
        print "Directory "+args.directory+" does not exist. I will try to create it"
        try:
            os.mkdir(args.directory)
        except:
            print "Failed to make directory "+args.directory
            return -1
    """
    
    if not os.path.isdir(args.directory):
        print "Directory "+args.directory+" does not exist."
        return -1
    
    if not os.path.isdir(args.directory+"/"+args.genus):
        os.mkdir(args.directory+"/"+args.genus)
        
    
    ref_seqs = SeqIO.parse(args.reference,'fasta')
    
    out_f_handles = {}
    
    for sr in ref_seqs:
        m = header_re.search(sr.description)
        if m:
            taxonomy = m.group('taxonomy').split(';')
            genus = taxonomy[-1]
            organism = " ".join(m.group('organism').split(' ')[:2]) # Get rid of strain junk
            
            if args.genus == genus:
                print organism
                
                if not args.mock:
                    try:
                        of = out_f_handles[organism]
                    except KeyError:
                        of = open(args.directory+'/'+args.genus+'/'+organism+".fasta",'a')
                        out_f_handles[organism] = of
                        
                    SeqIO.write(sr,of,'fasta')
                else:
                    sys.stdout.write(sr.description+'\n')

    for org in out_f_handles:
        out_f_handles[org].close()
                
            
            
if __name__ == "__main__":
    main()