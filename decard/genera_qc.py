#!/usr/bin/env python
import argparse

import os
from Bio import SeqIO
import tempfile
from deenurp.subcommands import filter_outliers


def no_ambiguous_bases(seq):
    unambiguousBases = ['a','c','t','g','u']
    return (  len([b for b in seq.lower() if not b in unambiguousBases]) == 0)

def main():
    args_parser = argparse.ArgumentParser()
    args_parser.add_argument('genus_dir')
    args_parser.add_argument('--minimum_length','-m',default=1000)
    args_parser.add_argument('--percentile','-p',default=90)
    args_parser.add_argument('--check_min_n','-cm',default=3)
    
    args = args_parser.parse_args()
    
    minimum_length = int(args.minimum_length)
    
    # 1. Get all the files in this genus's directory, filtering for FASTA files. 
    genus_fasta = [fn for fn in os.listdir(args.genus_dir) if '.fasta' in fn.lower()]
    
    """
        2. We should make an effort to filter out dubiously named files here. Things that do not start with the genus name and look like crap
            Also need to treat the grab-bag Genus sp or Genus spp file that is likely made up of a mix of novel species. 
    """
    # 3. Open each file
    for fn in genus_fasta:
        print fn
        sp_sr = SeqIO.parse(args.genus_dir+"/"+fn,'fasta')
        
        # An associative dict for sequences that past muster
        species_sr = {}
        # 4. And loop through it's records
        for sr in sp_sr:
            # i. Check to see if the sequence length is over some goal minimum
            if len(sr.seq) < minimum_length:
                print sr.id, "Does not meet minimum length at ", len(sr.seq)
                continue
            # ii. Check for ambiguous bases. Don't bother if there are
            if not no_ambiguous_bases(sr.seq):
                print sr.id, "Has ambiguous bases"
                continue
            
            
            # This assignment also takes care of multiple sequences with the same ID
            # The last sequence of sufficient quality of this ID in the file is the one that is kept,
            species_sr[sr.id] = sr
        
        # 5. See how many sequences made the cut at this point. If less than one, we can delete the fasta file and move on. 
        if len(species_sr) < 1:
            print "Removing", fn
            os.remove(args.genus_dir+"/"+fn)
            continue

            
        # 6. Use the deenurp package to cluster and determine if there are outlier sequences         
        if len(species_sr) >= int(args.check_min_n):
        
            with tempfile.NamedTemporaryFile() as filtered_tf:
                for seq_id in species_sr:
                    SeqIO.write(species_sr[seq_id],filtered_tf,'fasta')
                filtered_tf.flush()
                
                filtered_df = filter_outliers.filter_sequences('',
                                                 sequence_file=filtered_tf.name,
                                                taxa=fn,
                                                strategy='cluster',
                                                #percentile=args.percentile,
                                                percentile=90,
                                                min_radius=0.0,
                                                max_radius=None,
                                                cluster_type='single',
                                                aligner='cmalign',
                                                executable=None,
                                                 )
                filtered_tf.close()
            

            not_outliers = filtered_df[filtered_df.is_out==False].seqname

        else:
            not_outliers = species_sr.keys()
            
        # Again, if nobody made the cut, delete the file
        if len(not_outliers) < 1:
            print "Removing", fn
            os.remove(args.genus_dir+"/"+fn)
            continue
        # Now save back to the file, removing what was there before
        with open(args.genus_dir+"/"+fn,'w') as species_f:
            for seq_id in not_outliers:
                SeqIO.write(species_sr[seq_id],species_f,'fasta')
            species_f.close()


if __name__ == "__main__":
    main()