#!/usr/bin/env python

import argparse
from Bio import SeqIO


def main():
    args_parser = argparse.ArgumentParser()

    args_parser.add_argument('fastq_in',  help='Input: sequence file in fastq format')
    args_parser.add_argument('fastq_out',  help='OUTPUT: Where to put the sequences, in fastq format (with stripped IDs)')
    
    args = args_parser.parse_args()
    
    try:
        fastq_recs = SeqIO.parse(args.fastq_in,'fastq')
    except:
        print "Could not open fastq file for reading"
        return -1

    try:
        fastq_out_h = open(args.fastq_out,'w')
    except IOError:
        print "Could not open fastq file for writing"
        return -1
    
    
    for fastq_r in fastq_recs:
        fastq_r.id = fastq_r.id[:-4]
        fastq_r.description = fastq_r.id
        fastq_r.name = fastq_r.id
        SeqIO.write(fastq_r,fastq_out_h,'fastq')
    
    fastq_out_h.close()
    
            
if __name__ == "__main__":
    main()
