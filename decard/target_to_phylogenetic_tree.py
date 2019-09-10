#!/usr/bin/env python
"""
        Goal: Take the target input, and covert it to something that can be dragged through the R
        Phyloseq package.
        Input: Targets.csv
        
"""
import argparse
import pandas as pd
import re
from Bio import SeqIO, AlignIO
from Bio import Entrez
import tempfile
import sys
import subprocess
import os

header_re = re.compile('description=\"(?P<description>[^\"]+)\" organism=\"(?P<organism>[^\"]+)\" taxonomy=\"(?P<taxonomy>[^\"]+)\" ncbi_tax_id=\"(?P<ncbi_tax_id>[^\"]+)\"')

def get_taxonomy_from_source(row):
        try:
                seq_recs = SeqIO.parse(row.source_file,'fasta')
        
                for sr in seq_recs:
                    if sr.id == row.sequence_id:
                        m = header_re.search(sr.description)
                        if m:
                            return m.group('taxonomy')
                        else:
                            return ""
                        break
                if not 'ncbi_tax_id' in target:
                    return None
        except IOError:
                return ""
def lookup_rank(taxonomy_list, rank_index):
        try:
                return taxonomy_list[rank_index]
        except IndexError:
                return None

def main():
        args_parser = argparse.ArgumentParser()
        
        args_parser.add_argument('Target', help='INPUT  (target) file, in CSV format')
        args_parser.add_argument('--output_dir', '-o', help='Where to put the output tree (newick)', required=True)
        args_parser.add_argument('--email', '-e', help="Email for use with NCBI", required = True)
        args_parser.add_argument('--cm', '-c', help='Covariance matrix for cmalign', required=True)
        args_parser.add_argument('--threads', '-t', help='Number of threads', default = 1)
        
        
        args = args_parser.parse_args()
        
        threads = str(int(args.threads))
        cm_file = args.cm
        
        output_dir = args.output_dir
        
        # Load data
        target_df = pd.read_csv(args.Target)
        
        temp_fasta_fh = tempfile.NamedTemporaryFile(mode='w', delete=False)
        
        # Some accessions may be duplicated. Use a set to transfer them only once. 
        transferred_accessions = set()
        # If our file is missing, put it here so we can grab them
        missing_accessions = set()
        # For each row, get the file and then sequence, put in into the temp fasta file
        for idx, row in target_df.iterrows():
                try:
                        source_file_reads = SeqIO.parse(row.source_file,'fasta')
                        for sr in source_file_reads:
                                if sr.id == row.sequence_id:
                                        if not sr.id in transferred_accessions:
                                                print "Transferring ", sr.id
                                                SeqIO.write(sr,temp_fasta_fh,'fasta')
                                                transferred_accessions.add(sr.id)
                                        # Now we can stop our loop and go to the next row  
                                        break
                except IOError:
                        sys.stderr.write("Could not find "+str(row.source_file)+"\n")
                        missing_accessions.add(row.sequence_id)
                        pass # No sense in searching in an empty file
        
        # OK, now try to grab the missing accessions from the NCBI nt database.
        
        for missing_acc in list(missing_accessions):
                
                try:
                        h = Entrez.efetch(db='nucleotide', id = missing_acc, retmode='fasta', rettype='fasta')
                except:
                        # accession is frequently in the format accc.start.end. Try that now
                        acc_list = missing_acc.split('.')
                        try:
                                start = int(acc_list[1])
                        except IndexError:
                                start = None
                        try:
                                end = int(acc_list[2])
                        except IndexError:
                                end = None
                        
                        try:
                                h = Entrez.efetch(db='nucleotide', id = acc_list[0], retmode='fasta', rettype='fasta')
                        except:
                                # Still no dice
                                h = None
                        if h:
                                try:
                                        sr = SeqIO.read(h,'fasta')
                                        # If we have a start and end, slice the recovered sequence to them
                                        if start and end:
                                                sr.seq = sr.seq[start:end]
                                        missing_accessions.remove(missing_acc)
                                        SeqIO.write(sr, temp_fasta_fh,'fasta')
                                except:
                                        pass
                                
        if len(missing_accessions) > 0:
                sys.stderr.write("Still missing accessions "+", ".join(list(missing_accessions))+"\n")
                
        # Be sure we flush everything to disk
        temp_fasta_fh.close()
        # Great. Now we need to align these sequences into phylip format. We will use cmalign
        
        # RUN
        # cmalign --cpu=4 --dnaout --noprob -o temp.sto /fh/fast/fredricks_d/pipeline_group/communities/data/RRNA_16S_BACTERIA.cm seqs.fasta
        
        
        aln_temp_name = next(tempfile._get_candidate_names())
        aln_temp_path = os.path.join(tempfile.gettempdir(), aln_temp_name)
        
        subprocess.call(['cmalign',
                         '--cpu='+threads,
                         '--dnaout',
                         '--noprob',
                         '-o'+aln_temp_path,
                         cm_file,
                         temp_fasta_fh.name,  
                        ])
        
        # Output is in stockholm. Convert to 'phylip-relaxed' with alignIO
        aln_phyl_h = tempfile.NamedTemporaryFile(mode='w',suffix='.phl', delete=False)
        aln = AlignIO.read(aln_temp_path,'stockholm')
        AlignIO.write(aln, aln_phyl_h, 'phylip-relaxed')
        aln_phyl_h.close()
        
        # Now make a tree with the alignment
        
        # raxmlHPC-PTHREADS-SSE3 -m GTRGAMMA -p 12345 -s align.phylip -n source.tre -T 4
        
        subprocess.call(['raxmlHPC',
                         '-mGTRGAMMA',
                         '-p12345',
                         '-T'+threads,
                         '-s'+aln_phyl_h.name,
                         '-w'+os.path.abspath(output_dir),
                         '-ntarget.tre'
                         ])
                        
                        
                        
                
        

if __name__ == '__main__':
    main()
