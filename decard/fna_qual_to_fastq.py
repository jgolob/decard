#!/usr/bin/env python

import argparse
from Bio import SeqIO, SeqRecord
import gzip


def main():
    args_parser = argparse.ArgumentParser()

    args_parser.add_argument('fna',  help='Input: FNA')
    args_parser.add_argument('qual', help='Input: qual')
    args_parser.add_argument('fastq',  help='OUTPUT: FastQ filename')
    args_parser.add_argument('--gzip','-g', action='store_true',help='Will gzip compress output')
    
    args = args_parser.parse_args()
    
    if args.gzip:
        out_h = gzip.open(args.fastq,'w')
    else:
        out_h = open(args.fastq,'w')
    
    for s, q in zip(SeqIO.parse(args.fna,'fasta'), SeqIO.parse(args.qual,'qual')):
        if s.id == q.id and len(s.seq) == len(q.letter_annotations['phred_quality']):
            sr = SeqRecord.SeqRecord(s.seq,id=s.id, description=s.description, letter_annotations=q.letter_annotations)
            SeqIO.write(sr,out_h,'fastq')
        else:
            print "Mismatched ", sr, qual
    
    out_h.close()

            
if __name__ == "__main__":
    main()
