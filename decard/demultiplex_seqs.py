#!/usr/bin/env python

import os
import pandas as pd
from Bio import SeqIO

def main():
    args_parser = argparse.ArgumentParser()

    args_parser.add_argument('seq_info', help='')
    args_parser.add_argument('seqs', help='')
    args_parser.add_argument('outdir', help='')
    args_parser.add_argument('--suffix', '-s')
    args_parser.add_argument('--prefix', '-p')
    
    seq_info_df = pd.read_csv(args.seq_info, header=-1, names=['seq', 'community'])
    
    ourdir = ars.outdir
    if args.suffix:
        suffix = args.suffix
    else:
        suffix = ""
    
    if args.prefix:
        prefix = args.prefix
    else:
        prefix = ""
    
    out_f_h = {}
    for community in seq_info_df.community.unique():
        out_f_h[community]= open(os.path.join(outdir,prefix+community+suffix+".fasta"),'w')
    
    seq_info_df.set_index('seq', inplace=True)
    
    seqs = SeqIO.parse(args.seqs,'fasta')
    for seq in seqs:
        SeqIO.write(seq,out_f_h[seq_info_df.loc[seq.id].community],'fasta')
        
    for out_h in out_f_h:
        out_f_h[out_h].close()
        
    
    
    
    
    
    
    
    
    
    
if __name__ == '__main__':
    main()