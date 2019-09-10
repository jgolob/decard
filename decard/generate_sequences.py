#!/usr/bin/env python
"""
This module is used to take a target file generated by generate_targets, and use it to run in-silico PCR to generate sequences 
"""



import argparse
from Bio import SeqIO, pairwise2, Seq
from Bio.SeqRecord import SeqRecord
import re
import os
import csv
import numpy as np
import random
import uuid


header_re = re.compile('description=\"(?P<description>[^\"]+)\" organism=\"(?P<organism>[^\"]+)\" taxonomy=\"(?P<taxonomy>[^\"]+)\" ncbi_tax_id=\"(?P<ncbi_tax_id>[^\"]+)\"')


def read_targets(targets_fn):
    targets = []
    with open(targets_fn) as targets_f:
        targets_reader = csv.DictReader(targets_f)
        
        targets = [t for t in targets_reader]
        
        targets_f.close()
    return targets


def read_primers(fasta_fn):
    seqs = SeqIO.parse(fasta_fn, "fasta")
    forward_primers = []
    reverse_primers = []
    for sr in seqs:
        if "Forward" in sr.id:
            forward_primers.append(sr.seq)
        elif "Reverse" in sr.id:
            reverse_primers.append(sr.seq)
    
    return forward_primers, reverse_primers

def in_silico_pcr(primers, seq, min_percent = 80, min_score=None):
    """
        A quick and dirty in-silico PCR method, using the pairwise alignment package from BioPython
        For now, it uses a percent cutoff. Ideally, it would use Tm-based modeling. (tbd....)
        
        Inputs:
            primers: A LIST of primers
            seq:    Sequence to try pcr against
            min_percent: Minimum percent match to consider this a match
            min_score: A minimum alignment score to consider acceptable.
                By default, the scoring is +1 for every match, -1 for every gap
    """

    hits = []
    for primer in primers:
        primer_len = len(primer)
        alns = pairwise2.align.localxs(primer,seq, -1, -1) + pairwise2.align.localxs(primer.reverse_complement(),seq, -1,-1)
        for aln in alns:
            match = str(aln[0]).replace("-","")
            score = aln[2]
            begin = aln[3]
            end  = aln[4]
            
            # First check to see if the match end corresponds to the end of the primer
            # If there isn't a 3' match, no need to continue
            if (str(primer)[-1:] == match[-1:]) or (str(primer.reverse_complement())[-1:] == match[-1:]):
                if min_percent and min_score:
                    if (score / float(primer_len) >= min_percent/100.0) and (score >= min_score):
                        hits.append([begin,end])
                elif min_percent:
                    if (score / float(primer_len) >=min_percent/100.0):
                        hits.append([begin,end])
                elif min_score:
                    if (score >= min_score):
                        hits.append([begin,end])
                else:
                    hits.append([begin,end])

    
    if len(hits) > 1:
        start = min([h[0] for h in hits])
        end = max([h[1] for h in hits])
        
        return seq[start:end]
    else:
        return None
                    

def generate_seq_id(prefix=None):
    if prefix:
        seq_id = prefix
    else:
        seq_id=""
    
    return seq_id+str(uuid.uuid4()).replace('-','')
    

def main():
    args_parser = argparse.ArgumentParser()

    args_parser.add_argument('--targets', '-t', help='CSV file with the targets for this community', required=True)
    args_parser.add_argument('--primers','-p', help="Fasta file with primers", required=True)
    args_parser.add_argument('--rRNA16SSU','-r', help="Directory containing 16SSU rRNA", required=True)
    args_parser.add_argument('--output','-o', help="Fasta file to put our simulated reads into", required=True)
    args_parser.add_argument('--count', '-c', help='Total count of sequences per community', default = 5000)
    args_parser.add_argument('--map','-M', help="CSV file to map our reads to communities", required=True)
    args_parser.add_argument('--mock', '-m', help='Mock run. Do not modify the FASTA file and limit how many records we go after', action='store_true')
    
    
    args = args_parser.parse_args()
    count = int(args.count)
    
    rRNA16SSU = args.rRNA16SSU
    
    
    # Get our files ready. 
    if args.mock:
        import sys
        out_f = sys.stdout
        out_map_f = sys.stdout
    else:
        out_f = open(args.output,'w')
        out_map_f = open(args.map,'w')
    
    map_writer = csv.DictWriter(out_map_f, ['seqID','community','organism','ncbi_tax_id','sourceSeq'] )
    map_writer.writeheader()
        
    targets = read_targets(args.targets)
    
    forward_primers, reverse_primers = read_primers(args.primers)
    primers = forward_primers+reverse_primers

    for target in targets:
        
        reads_for_target = int(float(target['weight'])*count)
        if reads_for_target > 0:
            try:
                srs = SeqIO.parse(os.path.join(rRNA16SSU,target['source_file']),'fasta')
                found = False
                for sr in srs:
                    if sr.id == target['sequence_id']:
                        found = True
                        break
            except:
                print "Could not get ", target['sequence_id'],'from ', target['source_file']
                continue
            if found == False:
                continue
            
            seq = sr.seq
            seq = Seq.Seq(str(seq).replace('U','T'))
            # Run in-silico PCR
            amplicon = in_silico_pcr(primers, seq)
    
            # Maybe our primers won't amplify. 
            if not amplicon:
                print "Could not amplify from ", target['source_file']
            else:
            # Get the taxonomy ID
                m = header_re.search(sr.description)
                if m:
                    organism = m.group('organism')
                    ncbi_tax_id = m.group('ncbi_tax_id')
                else:
                    organism = target['organism']
                    ncbi_tax_id = None
                for n in xrange(reads_for_target):
                    curID = generate_seq_id(prefix=target['community_name']+'SCR')
                    amplicon_sr = SeqRecord(amplicon, id=curID, description='From '+str(sr.description))
                    SeqIO.write(amplicon_sr,out_f,'fasta')
                        
                    map_writer.writerow({
                        'community': target['community_name'],
                        'seqID':      curID,
                        'organism':     organism,
                        'ncbi_tax_id':  ncbi_tax_id,
                        'sourceSeq':    sr.id,
                        })

                        
        
        
    out_f.close()
    out_map_f.close()
               
            
if __name__ == "__main__":
    main()
