#!/usr/bin/env python

import argparse
from Bio import SeqIO
import re

re_description = re.compile('^(?P<badid>\S+)\s+(?P<oldid>\S+)\s+')

def fix_id(record):
    match = re_description.match(record.description)
    if match:
        return match.group('oldid')
    else:
        print "Could not find the old id in ", sr
        return record.id

def main():
    args_parser = argparse.ArgumentParser()

    args_parser.add_argument('pre',  help='Input: FNA after preprocessing in QIIME')
    args_parser.add_argument('post',  help='OUTPUT: FNA to allow by-sequence annotation')
    
    args = args_parser.parse_args()
    
    try:
        pre_recs = SeqIO.parse(args.pre,'fasta')
    except:
        print "Could not open fasta file for reading"
        return -1

    try:
        post_h = open(args.post,'w')
    except IOError:
        print "Could not open fasta file for writing"
        return -1
    
    for rec in pre_recs:
        
        rec.id = fix_id(rec)
        SeqIO.write(rec,post_h,'fasta')

    
    post_h.close()

            
if __name__ == "__main__":
    main()
