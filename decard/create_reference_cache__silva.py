#!/usr/bin/env python

import argparse
from Bio import SeqIO
import re
import os
import sys
import sqlite3
import gzip
from TaxonCache import TaxonCache

re_silva_description = re.compile('(?P<accession>[\w]+)\.(?P<acc_loc>[\d\.]+) (?P<taxonomy>[\w; ]+)')


def main():
    args_parser = argparse.ArgumentParser()
    args_parser.add_argument('--reference','-r', help='FASTA reference file', required=True)
    args_parser.add_argument('--database','-d', help='SQLite3 DB file', required=True)
    args_parser.add_argument('--mock', '-m', help='Mock run. Do not modify the FASTA file and limit how many records we go after', action='store_true')
    args_parser.add_argument('--email','-e',required=True)
    args_parser.add_argument('--url','-u',help='Base URL to REST db with taxonomic information. If not provided, we will default to NCBI')
    
    args = args_parser.parse_args()
    
    if args.url:
        taxons = TaxonCache(source='custom', url = args.url, email=args.email)
    else:
        taxons = TaxonCache(email=args.email)
        
    try:
        ref_f = gzip.open(args.reference)
        ref_seqs = SeqIO.parse(ref_f,'fasta')
    except:
        print "Could not open reference fasta file"
        return -1
    
    try:
        reference_cache = sqlite3.connect(args.database)
        reference_cache_cursor = reference_cache.cursor()
    except:
        print "Could not open sqlite3 database (or create)"
        return -1
    
    # Prep the database
    try:
        reference_cache_cursor.execute("""CREATE TABLE accession_calls
                                   (accession_label PRIMARY KEY,
                                   organism,
                                   ncbi_tax_id,
                                   taxonomy,
                                   silva_taxonomy)""")
    except sqlite3.OperationalError:
        pass # Table probably already exists
    
    try:
        reference_cache_cursor.execute("""CREATE TABLE ncbi_taxon
                                   (ncbi_tax_id PRIMARY KEY,
                                   scientific_name,
                                   'rank',
                                   'parent_ncbi_tax_id',
                                   'color',
                                   'ncbi_taxonomy'
                                   )""")
    except sqlite3.OperationalError:
        pass # Table probably already exists

    try:
        reference_cache_cursor.execute("""CREATE TABLE silva_to_ncbi_taxon (ncbi_tax_id, silva_taxonomy_string)""")
    except sqlite3.OperationalError:
        pass # Table probably already exists
    

    for sr in ref_seqs:
        m = re_silva_description.search(sr.description)
        if m:
            # First try to look up in our DB
            taxonomyString = m.group('taxonomy')
            reference_cache_cursor.execute("""SELECT tax.scientific_name, tax.ncbi_tax_id, tax.rank, tax.ncbi_taxonomy
                                           FROM silva_to_ncbi_taxon, ncbi_taxon as tax
                                           WHERE silva_to_ncbi_taxon.silva_taxonomy_string = ? AND tax.ncbi_tax_id=silva_to_ncbi_taxon.ncbi_tax_id""", (taxonomyString,))
            cached_tax_match = reference_cache_cursor.fetchone()
            if cached_tax_match:
                organism, ncbi_tax_id, rank, ncbi_taxonomy = cached_tax_match
                
            else: # Don't have this taxon yet
                taxonomyList = taxonomyString.split(';')
                tax = taxons.lookup_taxonomy(taxonomyList)
                if tax:
                    ncbi_tax_id = tax["TaxId"]
                    organism = tax['ScientificName']
                    rank    = tax['Rank']
                    parent_tax_id = tax['ParentTaxId']
                    ncbi_taxonomy = '; '.join([t["ScientificName"] for t in tax['LineageEx']])
                    try:
                        reference_cache_cursor.execute("""INSERT INTO ncbi_taxon
                                                       (ncbi_tax_id, scientific_name, rank, parent_ncbi_tax_id, ncbi_taxonomy)
                                                       VALUES (?,?,?,?,?)
                                                       """, (ncbi_tax_id, organism, rank, parent_tax_id, ncbi_taxonomy))
                    except sqlite3.IntegrityError:
                        pass
                    try:
                        reference_cache_cursor.execute("""INSERT INTO silva_to_ncbi_taxon
                                                       (ncbi_tax_id, silva_taxonomy_string)
                                                       VALUES (?,?)
                                                       """, (ncbi_tax_id, taxonomyString))
                    except sqlite3.IntegrityError:
                        pass
                else:
                    continue
                    
                    
    
        
            try:
                reference_cache_cursor.execute("""INSERT INTO accession_calls
                                                    (organism, ncbi_tax_id, taxonomy, accession_label,silva_taxonomy)
                                                    VALUES
                                                    (?,?,?,?,?)
                                                """, (organism, ncbi_tax_id, ncbi_taxonomy, sr.id, taxonomyString) )
                
            except sqlite3.IntegrityError: # Entry already in the database, update instead
                reference_cache_cursor.execute("""UPDATE accession_calls SET
                                                    organism=?,
                                                    ncbi_tax_id=?,
                                                    taxonomy=?,
                                                    silva_taxonomy=?
                                                    WHERE accession_label = ?
                                                """, (organism, ncbi_tax_id, ncbi_taxonomy, taxonomyString, sr.id) )
            reference_cache.commit()
            
                
            
            
if __name__ == "__main__":
    main()