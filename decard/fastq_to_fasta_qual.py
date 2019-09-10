#!/usr/bin/env python

import argparse
from Bio import SeqIO


def main():
    args_parser = argparse.ArgumentParser()

    args_parser.add_argument('fastq',  help='Input: sequence file in fastq format')
    args_parser.add_argument('fasta',  help='OUTPUT: Where to put the sequences, in fasta format')
    args_parser.add_argument('qual',  help='OUTPUT: Where to put the quality scores, in qual format')
    
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
    
    try:
        qual_h = open(args.qual,'w')
    except IOError:
        print "Could not open qual file for writing"
        return -1
    
    for fastq_r in fastq_recs:
        SeqIO.write(fastq_r,fasta_h,'fasta')
        SeqIO.write(fastq_r,qual_h,'qual')
    
    qual_h.close()
    fasta_h.close()
    
            
if __name__ == "__main__":
    main()
