#!/usr/bin/env python

import argparse
from Bio import SeqIO


def main():
    args_parser = argparse.ArgumentParser()

    args_parser.add_argument('fastq',  help='Input: sequence file in fastq format')
    args_parser.add_argument('fasta',  help='OUTPUT: Where to put the sequences, in fasta format')
    args_parser.add_argument('--no_ambiguous','-na', help='Option to drop sequences with ambiguous', action='store_true')
    
    args = args_parser.parse_args()
    
    try:
        fastq_recs = SeqIO.parse(args.fastq,'fastq')
    except:
        print "Could not open fastq file for reading"
        return -1

    try:
        fasta_h = open(args.fasta,'w')
    except IOError:
        print "Could not open fasta file for writing"
        return -1
    
    
    for fastq_r in fastq_recs:
        if args.no_ambiguous:
            if len([s for s in str(fastq_r.seq) if s not in ['A','a','G','g','C','c','T','t','U','u']]) == 0:
                SeqIO.write(fastq_r,fasta_h,'fasta')
        else:
            SeqIO.write(fastq_r,fasta_h,'fasta')
        
    fasta_h.close()
    
            
if __name__ == "__main__":
    main()
