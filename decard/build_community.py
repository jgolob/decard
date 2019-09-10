#!/usr/bin/env python

import argparse
from Bio import SeqIO, pairwise2
from Bio.SeqRecord import SeqRecord
import re
import os
import csv
import numpy as np
import random
import uuid

header_re = re.compile('description=\"(?P<description>[^\"]+)\" organism=\"(?P<organism>[^\"]+)\" taxonomy=\"(?P<taxonomy>[^\"]+)\" ncbi_tax_id=\"(?P<ncbi_tax_id>[^\"]+)\"')

def read_distro_goal(distribution_fn):
    """
        Given a filename, read in the distribution goal from a properly formatted CSV
        Should be:
        Genus, Fraction, STD, Species_n, Species_std, Species_slope, Species_intercept, Species_p
        
    """
    goal = []
    with open(distribution_fn) as distro_f:
        reader = csv.DictReader(distro_f)
        for row in reader:
            for k in row:
                try:
                    row[k]=float(row[k])
                except:
                    pass
            goal.append(row)
        distro_f.close()
        
    return goal

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

def in_silico_pcr(primers, seq, min_percent = 90, min_score=None):
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
            score = aln[2]
            begin = aln[3]
            end  = aln[4]
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
    
    if len(hits) > 0:
        start = min([h[0] for h in hits])
        end = max([h[1] for h in hits])
        
        return seq[start:end]
    else:
        return None
                    
    
def communities_generateGoals(num, distribution_goal, count):
    communities = []
    
    for i in range(0,int(num)):
        community = {
            'num':  i,
            'name': "SYNCM"+str(i)+"_"+str(uuid.uuid4()).replace('-','')
        }
        print "_____ ", i, " ________"
        community_fract = []
        # For each goal fraction and std, generate a random value based on the normal distribution for what proportion the sample should be for this community
        for genus_goal in distribution_goal:
            # Use the typical fractional abundance of this genus, plus the std of this mean to give us a random fraction (abundance)
            fract = np.max([0.0,np.random.normal(loc = genus_goal['Fraction'], scale = genus_goal['STD'])])
            
            if fract > 0.0:
                # Then decide how many species to pull (richness)
                
                # If the log regression is good enough, use it
                
                if float(genus_goal['Species_p'] <= 0.05):
                    species = int(np.round(np.max([1, float(genus_goal['Species_slope']*np.log(fract)+float(genus_goal['Species_intercept']))])))
                # If the log model isn't great, see if we have a std deviation to use with the mean number to get our next estimate.
                elif int(genus_goal['Species_std']) > 0:
                    species = int(np.round(np.max([1, np.random.normal(loc=genus_goal['Species_n'], scale=genus_goal['Species_std'])])))
                else: # Just use the mean number
                    species = np.max([int(genus_goal['Species_n']),1])
                
                community_fract.append({
                                'genus': genus_goal['Genus'],
                                'fraction': fract,
                                'species_n': species})
        # Get the total proportion to be able to normalize
        total_fract = np.sum([cf['fraction'] for cf in community_fract])

        # Figure out the goal number of sequences per genera, using the proportions (normalized) and the goal total number of sequenesc
        genus_seqs = [{'genus': cf['genus'], 'count': int(np.round(cf['fraction']/total_fract*(count))), 'species_n': cf['species_n'] } for cf in community_fract ]
        genus_seqs.sort(key=lambda g: -g['count'])
        for gs in genus_seqs:
            print gs['genus'], gs['count'], gs['species_n']
        community['goal']= genus_seqs
        communities.append(community)

        
    return communities


def generate_seq_id(prefix=None):
    if prefix:
        seq_id = prefix
    else:
        seq_id=""
    
    return seq_id+str(uuid.uuid4()).replace('-','')
    

def main():
    args_parser = argparse.ArgumentParser()

    args_parser.add_argument('--genera_fasta', '-g', help='Directory where to find genera', required=True)
    args_parser.add_argument('--distribution','-d', help='CSV file with desired distribution, by genus', required=True)
    args_parser.add_argument('--mock', '-m', help='Mock run. Do not modify the FASTA file and limit how many records we go after', action='store_true')
    args_parser.add_argument('--number', '-n', help='How many communities to generate', default = 1)
    args_parser.add_argument('--count', '-c', help='Total count of sequences per community', default = 5000)
    args_parser.add_argument('--primers','-p', help="Fasta file with primers", required=True)
    args_parser.add_argument('--output','-o', help="Fasta file to put our simulated reads into", required=True)
    args_parser.add_argument('--map','-M', help="CSV file to map our reads to communities", required=True)
    
    args = args_parser.parse_args()
    
    count = int(args.count)
    
    if not args.mock:
        out_f = open(args.output,'w')
        out_map_f = open(args.map,'w')
        map_writer = csv.DictWriter(out_map_f, ['seqID','community','organism','ncbi_tax_id','sourceSeq'] )
        map_writer.writeheader()
        
              
    if not os.path.isdir(args.genera_fasta):
        print "Directory for genus fasta files "+args.genera_fasta+" does not exist"
        return -1
    
    distribution_goal = read_distro_goal(args.distribution)
    
    forward_primers, reverse_primers = read_primers(args.primers)
    primers = forward_primers+reverse_primers
    
    communities = communities_generateGoals(args.number, distribution_goal, count)
    
    cached_PCR = {}
    genera_dir = args.genera_fasta
    
    if not args.mock:
    
        for comm in communities:
            for genus_goal in comm['goal']:
                print "Loading sequences for ", genus_goal['genus']
                if not os.path.isdir(genera_dir+"/"+genus_goal['genus']):
                    print "Missing "+genus_goal['genus']+"'s directory"
                else:
                    # Get the filenames of each species fasta file in this genus
                    fasta_fns =  os.listdir(genera_dir+"/"+genus_goal['genus'])
                    # Randomize the order of the filenames
                    random.shuffle(fasta_fns)
                    # Can't take more species than we have, so take the min of the two
                    species_n = min(len(fasta_fns), genus_goal['species_n'])
                    # Cut the list down to the wanted number of species
                    fasta_fns = fasta_fns[:species_n]
                    for i in range(genus_goal['count']):
                    # Try once for each of our goal count
                        fasta_fn = random.choice(fasta_fns)
                        # Graba fasta file for this species
                        genus_srs = []
                        seqs = SeqIO.parse(genera_dir+'/'+genus_goal['genus']+'/'+fasta_fn,'fasta')
                        for sr in seqs:
                                # parse and load the records
                            genus_srs.append(sr)
                        
                        # Pick a random seq    
                        sr = random.choice(genus_srs)
                        
                        # See if we have a cached result for this PCR, as the in-silico PCR is SLOW
                        try:
                            amplicon = cached_PCR[sr.id]
                        except KeyError:
                            seq = sr.seq
                            # Run in-silico PCR
                            amplicon = in_silico_pcr(primers, seq)
                            cached_PCR[sr.id]=amplicon
                        # Maybe our primers won't amplify. Prinr o
                        if not amplicon:
                            print "Could not amplify from ", fasta_fn
                        else:
                            curID = generate_seq_id(prefix='SCR'+str(comm['num']))
                            amplicon_sr = SeqRecord(amplicon, id=curID, description='From '+fasta_fn+" "+str(sr.id))
                            SeqIO.write(amplicon_sr,out_f,'fasta')
                            # Get the taxonomy ID
                            m = header_re.search(sr.description)
                            if m:
                                organism = m.group('organism')
                                ncbi_tax_id = m.group('ncbi_tax_id')
                            else:
                                organism = fasta_fn.split('.')[0]
                                ncbi_tax_id = None
                                
                            map_writer.writerow({
                                'community': comm['name'],
                                'seqID':      curID,
                                'organism':     organism,
                                'ncbi_tax_id':  ncbi_tax_id,
                                'sourceSeq':    sr.id,
                                })
        
                        
        
        
        out_f.close()
        out_map_f.close()
        
                        
                
                
            
            
if __name__ == "__main__":
    main()
