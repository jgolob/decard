#!/usr/bin/env python
"""
        Goal: Take the target input, and covert it to something that can be dragged through the R
        Phyloseq package.
        Input: Targets.csv
        
"""
import argparse
import pandas as pd

#import json
#import re
from Bio import SeqIO
import sys
import math
import os


def main():
        args_parser = argparse.ArgumentParser()
        
        args_parser.add_argument('fasta', help='INPUT concatenated fasta file')
        args_parser.add_argument('output_prefix', help='Prefix (including path) to use for output fasta')
        args_parser.add_argument('--seq_info', '-si', help='Seq Info file mapping seq_id to community (headerless)', required=True)
        args_parser.add_argument('--group_size', '-g', help='OPTIONAL how many communities to put in each chunk. If negative, each community gets its own')
       
        
        args = args_parser.parse_args()
        
        
        seq_info_df = pd.read_csv(args.seq_info, header=None, names=['seq_id','community'])
        seq_info_df.index = seq_info_df.seq_id
        
        output_prefix = args.output_prefix
        
        communities = list(seq_info_df.community.unique())
        
        communities.sort()
        
        if args.group_size:
                group_size = int(args.group_size)
        else:
                group_size = 1
        
        
        num_chunks = int(math.ceil(len(communities) / float(group_size)))
        
        # A dict mapping a community to a particular chunk file handle
        community_to_chunk_fh = {}        
        
        for i in xrange(num_chunks):
                chunk_fh = open(output_prefix+str(i)+".fasta",'w')
                communities_in_chunk = communities[i*group_size:(i+1)*group_size]
                print len(communities_in_chunk), communities_in_chunk
                for community in communities_in_chunk:
                      community_to_chunk_fh[community] = chunk_fh
                      
        for sr in SeqIO.parse(args.fasta,'fasta'):
                community = seq_info_df.loc[sr.id].community
                SeqIO.write(sr, community_to_chunk_fh[community],'fasta')
                      
        
                        
                
        

if __name__ == '__main__':
    main()
