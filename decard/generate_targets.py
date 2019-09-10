#!/usr/bin/env python
"""
    This module generates specific targets for desired distribution, and a given set of reference sequences
    INPUTS:
        # distribution: In CSV format.
            Genus, Fraction, STD, Species_n, Species_std, Species_slope, Species_intercept, Species_p
        # Genera Directory:
            A directory in which each genus has it's own subdirectory.
            Within that directory is a fasta file for each species within that genus
            Each fasta file should have one or more representitive sequences
            eg:
            /path/to/genera_dir/Streptococcus/Streptococcus pyogenes.fasta
        # Count:
            Number of targets per community
        # Number of communities:
            How many communities to generate, with the given count and distribution provided.
        
        # Prefix (opt): Prefix to add to the start of the community name
        # Suffix (opt): Suffix to add to the end of the community name
    
    OUTPUTs (for use in generate_amplicons):
        # targets_csv: A CSV file with the following columns to store the generated community targets.
        community_name, source_file, species, sequence_id, weight
    
"""
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
        
        Outputs a list of dicts
        
    """
    goal = []
    with open(distribution_fn,'rU') as distro_f:
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
   
def communities_generateGoals(num, distribution_goal, distribution_name="", prefix="", suffix="",offset=0):
    communities = []
    
    for i in xrange(offset,offset+num):
        community = {
            'num':  i,
            'name': prefix+"CM"+str(i)+suffix,
            'distribution_name':    distribution_name,
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
        genus_seqs = [{'genus': cf['genus'], 'fraction': cf['fraction']/total_fract, 'species_n': cf['species_n'] } for cf in community_fract ]
        
        # Get rid of any that end up with a zero count
        genus_seqs = [gs for gs in genus_seqs if gs['fraction'] > 0]
        # Sort for prettiness
        
        genus_seqs.sort(key=lambda g: -g['fraction'])
        # Print for niceness to our user
        for gs in genus_seqs:
            print gs['genus'], gs['fraction'], gs['species_n']
        community['goal']= genus_seqs
        communities.append(community)

        
    return communities

def calculate_species_count(nth_species, total_num_species , total_count):
    """
        Goal: Use a log model to get a count for the nth species, given a total num of species, and overall count
        nth_species: 1 to total_num (starts at 1)
        total_num_species: total num of species for total_Count
        total_count: total num for this genus
        
        In log10 space, we want to scale down to 1 - 10. Log10(1) = 0, Log10(10) = 1. We can then take the difference of nth_species - n-1th species scaled to get the count
        The net result at the end due to rounding may be off. That's fine. 
    """
    """ First scale to [1,10], and get where we are on this scaled value
        Why? We want our values to vary from 1-10 linearly
        s = mX + b.
        we want s = 1 when X = 0 
        so b = 1
        we want s = 10 when X = MAX
        so solving for m
        10 = m(MAX) + 1
        9 = m(MAX)
        m = 9 / MAX
    """
    m = 9 / float(total_num_species)
    b = 1.0
    s = m * nth_species + b
    s0 = m*(nth_species-1) + b
    
    # Next log transform our S
    s_l = np.log10(s)
    s0_l = np.log10(s0)
    
    return total_count*(s_l- s0_l)

def communities_pickSequences(communities, genera_dir):
    
    # The list into which we will put our targeted sequences
    sequences = []
    for comm in communities:
        for genus_goal in comm['goal']:
            print "Loading sequences for ", genus_goal['genus']
            if not os.path.isdir(genera_dir+"/"+genus_goal['genus']):
                print "Missing "+genus_goal['genus']+"'s directory"
            else:
                # Get the filenames of each species fasta file in this genus
                fasta_fns =  os.listdir(genera_dir+"/"+genus_goal['genus'])
                if len(fasta_fns) < 1:
                    print "No species available for ", genus_goal['genus']
                else:
                    # Randomize the order of the filenames
                    random.shuffle(fasta_fns)
                    # Can't take more species than we have, so take the min of the two
                    species_n = min(len(fasta_fns), genus_goal['species_n'])
                    # Cut the list down to the wanted number of species
                    fasta_fns = fasta_fns[:species_n]
                    
                    # Great, for each file representing a sequence.....
                    for n,fasta_fn in enumerate(fasta_fns):
                        # Figure out how many copies of this species we should have
                        species_count = calculate_species_count(n+1, species_n, genus_goal['fraction']*100)
                        # Grab the file.....
                        seqs = SeqIO.parse(genera_dir+'/'+genus_goal['genus']+'/'+fasta_fn,'fasta')
                        # Load all the strains / representatives into an array
                        species_srs = []
                        seqs = SeqIO.parse(genera_dir+'/'+genus_goal['genus']+'/'+fasta_fn,'fasta')
                        for sr in seqs:
                            # load the records
                            species_srs.append(sr)
                            # Pick a random sequence for this species
                        sr = random.choice(species_srs)
                        # Add it to our list
                        sequences.append({
                            'community_name': comm['name'],
                            'distribution_name': comm['distribution_name'],
                            'source_file':    os.path.abspath(genera_dir+'/'+genus_goal['genus']+'/'+fasta_fn),
                            'species': fasta_fn.replace('.fasta',''),
                            'sequence_id':  sr.id,
                            'weight':   species_count/100,
                        })
                        
    return sequences
                                         
def main():
    args_parser = argparse.ArgumentParser()

    args_parser.add_argument('--genera_fasta', '-g', help='Directory where to find genera', required=True)
    args_parser.add_argument('--distribution','-d', help='CSV file(s) with desired distribution(s), by genus', nargs='*', required=True)
    args_parser.add_argument('--mock', '-m', help='Mock run. Do not modify the FASTA file and limit how many records we go after', action='store_true')
    args_parser.add_argument('--number', '-n', help='How many communities to generate (per distribution)', nargs='*', default = [1])
    args_parser.add_argument('--output','-o', help="Output file for targets for PCR step / generate amplicons", required=True)
    args_parser.add_argument('--prefix','-p', help="Prefix to prepend to community names", default="")
    args_parser.add_argument('--suffix','-s', help="Suffix to append to community names", default="")
    
    
    args = args_parser.parse_args()
    
    # First handle our distros
    distribution_files = args.distribution
    distribution_goals = [read_distro_goal(distro) for distro in distribution_files]
    
    
    # Unpack and tidy up our num per distro
    numbers = [int(num) for num in args.number]

    if len(numbers) != len(distribution_goals):
        if len(numbers) == 1:
            numbers = numbers*len(distribution_goals)
        else:
            print "Please match number of distributions to number of communities per distribution"
            return -1
    
    
    
    
    # See if our genera dir exists (and could do some validation testing too if we wanted)
    if not os.path.isdir(args.genera_fasta):
        print "Directory for genus fasta files "+args.genera_fasta+" does not exist"
        return -1
    # Implicit else we're good to go
    genera_dir = args.genera_fasta
    
    # If we're not in mock mode, output to files
    if not args.mock:
        out_f = open(args.output,'w')
    else: # We are mock, output to stdout
        import sys
        out_f = sys.stdout
    
    # Set up writers.
    target_writer = csv.DictWriter(out_f, ['community_name', 'distribution_name', 'source_file','species', 'sequence_id', 'weight'])
    target_writer.writeheader()
    
    communities = []
    offset = 0
    
    for i, (distribution_goal, number) in enumerate(zip(distribution_goals,numbers)):
        distro_name = distribution_files[i]
        communities = communities+communities_generateGoals(number, distribution_goal, distribution_name=distro_name, prefix=args.prefix, suffix=args.suffix, offset=offset)
        offset+=number
    
    sequences = communities_pickSequences(communities, genera_dir)
    
    target_writer.writerows(sequences)    
    out_f.close()         
            
if __name__ == "__main__":
    main()
