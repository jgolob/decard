#!/usr/bin/env python
"""
    This program downloads and then updates / maintains a local FASTA file of sequences from the
    NCBI 16S Microbial database, for use in classifications
"""


"""
    For downstream, to extract key variables from the sr.description:
    description_re = re.compile('description=\"(?P<description>[^\"]+)\"')
    organism_re = re.compile('organism=\"(?P<organism>[^\"]+)\"')
    taxonomy_re = re.compile('taxonomy=\"(?P<taxonomy>[^\"]+)\"')
    ncbi_tax_id_re = re.compile('ncbi_tax_id=\"(?P<ncbi_tax_id>[^\"]+)\"')
    
    header_re = re.compile('description=\"(?P<description>[^\"]+)\" organism=\"(?P<organism>[^\"]+)\" taxonomy=\"(?P<taxonomy>[^\"]+)\" ncbi_tax_id=\"(?P<ncbi_tax_id>[^\"]+)\"')
    
    
"""
import argparse
from Bio import Entrez, SeqIO

## CONSTANTS
NCBI_16S_MICROBIAL_TERM = '33175[BioProject]'


def get_current_refs(f):
    """
        Read the existing FASTA file, and store the IDs in a list.
        Loads the sequences into an associative DB.
        This is a bit wasteful for memory.
    """
    cur_seqs = {}

    seqs = SeqIO.parse(f,'fasta')
    
    for seq in seqs:
        cur_seqs[seq.id] = seq
    
    return cur_seqs

def get_current_NCBI_16S():
    """
        Search the NCBI NT database for the Microbial 16S project (term written above)
        Return a list of the current IDS
    """
    RETMAX = 100000
    
    search_h = Entrez.esearch(db='nucleotide', term=NCBI_16S_MICROBIAL_TERM, retmax=RETMAX)
    search_r = Entrez.read(search_h)
    """
        Retmax is an annoying part of this function we have to handle. If there are more records than we're allowed to acquire at once
        We will have to handle batching.
    """
    idList = search_r['IdList']
    if int(search_r['Count']) > RETMAX:
        blocks = int(search_r['Count']) / RETMAX
        for i in range(1,blocks+1):
            search_h = Entrez.esearch(db='nucleotide', term=NCBI_16S_MICROBIAL_TERM, retmax=RETMAX, retstart=(RETMAX*i))
            search_r = Entrez.read(search_h)
            idList=idList+search_r['IdList']
        
    #print len(idList) == len(set(idList)), len(idList)== int(search_r['Count']), len(idList), int(search_r['Count'])

    return idList

def match_existing_to_current(cur_seqs, ncbi_ids):
    """
        Given a list of our local copy's current IDs and the NCBI NT DB IDs,
        figure out the IDs missing from the local copy that are still in NCBI
        and those that are removed from NCBI but present in our local copy
    """
    missing_ids = [gbID for gbID in ncbi_ids if not gbID in cur_seqs]
    removed_ids = [gbID for gbID in cur_seqs if not gbID in ncbi_ids]
    
    return missing_ids, removed_ids
    
def download_from_ncbi(gbID):
    """
        Does as the name implies, download a sequence record from the NT database with the given genbank ID
        Return the record
    """
    
    handle = Entrez.efetch(db="nucleotide", id=gbID, rettype="gbwithparts", retmode="text")
    
    sr = SeqIO.read(handle,'genbank')
    
    # Here we will modifiy the description to aid our deconvolution later to an organism
    description = sr.description
    try:
        annotations = sr.annotations
    except AttributeError:
        print sr.id, ' missing annotation '
        annotations= None
        
    organism = annotations['organism']
    taxonomy = annotations['taxonomy']
    taxonomy_string = ';'.join(map(str, taxonomy))
    
    # And then a bit of fun to extract the taxID from this record
    
    try:
        source_f = [f for f in sr.features if f.type=='source'][0]
        taxID = [xref for xref in source_f.qualifiers['db_xref'] if 'taxon:' in xref][0].split('taxon:')[1]
    except:
        taxID = None
    
    sr.description='description=\"'+description+'\" organism=\"'+organism+'\" taxonomy=\"'+taxonomy_string+'\" ncbi_tax_id=\"'+str(taxID)+'\"'
    
    return sr 

def main():
    args_parser = argparse.ArgumentParser()

    args_parser.add_argument('output', help='Location of existing FASTA file, to be updated.')
    args_parser.add_argument('--email','-e', help='REQUIRED: Your Email, to send to NCBI.', required=True)
    args_parser.add_argument('--mock', '-m', help='Mock run. Do not modify the FASTA file and limit how many records we go after', action='store_true')
    
    args = args_parser.parse_args()
    
    Entrez.email = args.email
    
    try:
        existing_f = open(args.output,'r')
        cur_seqs = get_current_refs(existing_f)
        existing_f.close()
    except IOError: # The DB file didn't exist. Assume we need to make it fresh, and download everything
        cur_seqs = {}
    
    if args.mock:
        ncbi_ids = get_current_NCBI_16S()[:5]
    else:
        ncbi_ids = get_current_NCBI_16S()
        
    
    missing_ids, removed_ids = match_existing_to_current(cur_seqs, ncbi_ids)
    
    if args.mock:
        for gbID in ncbi_ids:
            try:
                    sr = cur_seqs[gbID]
            except KeyError: # Not in our current
                sr = download_from_ncbi(gbID)
            
            print sr.id, sr.description
            
    else:    
        with open(args.output,'w') as output_f:
        # Wipe out and prepare to create it anew
            for gbID in ncbi_ids:
                # First see if we have this sequence in our existing local records
                try:
                    sr = cur_seqs[gbID]
                except KeyError: # Not in our current
                    sr = download_from_ncbi(gbID)
                
                SeqIO.write(sr,output_f,'fasta')
            
            output_f.close()
                
if __name__ == "__main__":
    main()