#!/usr/bin/env python
import argparse
import pandas as pd
import numpy as np
import csv
import re
import json
from TaxonCache import TaxonCache

TARGET_RANKS = [None, 'species','genus']

def handle_tax_ids(raw_tax_id):
        # Some classifiers can narrow things down to a handful of tax ids. So we need to catch that sort of output
        # Convert everyone to tuples
        
        # First see if we have a value (vs nan)
        if pd.isnull(raw_tax_id):
                return tuple([np.nan])
        # Implicit else
        try:
            return tuple([int(raw_tax_id)])
        except ValueError: # we don't have a simple integer
            return tuple([int(tid) for tid in raw_tax_id.split(',')])    

def calculate_ranks_off(lineage_1, lineage_2, off = 0):
    if not lineage_1 and not lineage_2:
        return off
    if not lineage_1:
        return off + len(lineage_2)
    
    if not lineage_2:
        return off + len(lineage_1)
    
    if lineage_1[-1]==lineage_2[-1]:
        return off
    if len(lineage_1) == len(lineage_2):
        return 2 + calculate_ranks_off(lineage_1[:-1],lineage_2[:-1],off)
    elif len(lineage_1) < len(lineage_2):
        return 1 + calculate_ranks_off(lineage_1, lineage_2[:-1], off)
    else:
        return 1 + calculate_ranks_off(lineage_2, lineage_1[:-1], off)

def trim_and_unpack_taxon(taxon, target_rank):
        
        """
        Take a taxon object, as we'd recieve from NCBI (or a private cache),
        trim it to a desired rank, and return unpacked variables useful for subsequent processing.
        Input:
                A TAXON as a set of nested dictionaries
                target_rank:
                        In NCBI format (lower case. Not guaranteed for all)
                        IF none, will not bother with targeting the rank.
                        
        Outputs:
                taxonomy / lineage as an ordered list of taxon dicts
                taxonomy_tax_ids = an ordered list of tax ids (ints) for this lineage
                taxonomy_string= Semi-delimited string. Includes taxon name.
                
        
        Strategies:
                When targeting a rank, there are a few possibilities of where we're starting.
                1) The taxon is already at this rank. 
                        Proceed, without further processing
                2) The taxon is below the target rank:
                        TRIM taxon lineage until we're at this rank (provided it exists)
                3) The taxon is ABOVE this rank:
                        Proceed without further processing.
                        
                The algorithm:
                        Test to see the taxon is at this rank. If so, proceed.
                                If not:
                                        Test to see if this rank is in the taxon's lineage
                                                If so, trim to this taxon
                                                If not, continue
                                                
                                        
        """
        
        # Get our taxonomy / lineage as a list
        taxonomy = list(taxon['LineageEx'])
        # Check to see if this taxon ends our list. If not, add it
        if len(taxonomy) == 0:
                taxonomy=[taxon]
        elif not taxonomy[-1]['TaxId'] == taxon['TaxId']:
                taxonomy.append(taxon)
        
        if target_rank:
                # Test to see if we're already at our target rank
                if taxonomy[-1]['Rank'] == target_rank:
                        pass
                elif target_rank in [t['Rank'] for t in taxonomy]: #We aren't. Next test to see if the target rank is in this lineage
                        # Yup, we have this rank in our taxonomy. Slice our taxonomy to this rank
                        taxonomy = taxonomy[:[t['Rank'] for t in taxonomy].index(target_rank)+1]
                else:
                        """ We are not at the target rank, AND the target rank is not within our taxonomy
                            Two possiblities:
                                1. We are above our target rank.
                                2. We are below the target rank, but the target rank is not in our taxonomy.
                                        The NCBI taxonomy is unfortunate in that #2 is a possibility. Alas
                                        We could do some checking here and trim CLOSE to this taxonomy
                                        It's less of an issue for species or species groups. 
                        """
                        pass
                        
        return (
                taxonomy,
                [t['TaxId'] for t in taxonomy],
                "; ".join([t['ScientificName'] for t in taxonomy])
                
        )
        
        

def test_match(taxons, ncbi_tax_id__true, ncbi_tax_id__test, organism, name, target_rank=None):
        """
        Try to figure out if the classification matches the true organism
        for a given target rank
        INPUTS:
                taxons: A TaxonCache instance
                ncbi_tax_id__true: True NCBI tax ID
                ncbi_tax_id__test: Test NCBI tax ID
                organism: True organism name
                name: classified Name
                target_rank: Rank at which to match. If none, try to match perfectly.
                Suggested choices:
                        'species'
                        'genus'

        
        OUTPUT: A dict of
        
        {
        'result': ,
        'ranks_off':,
        'lineage_true':,
        'lineage_classified':,
        'taxon_true':,
        'taxon_classified':}
                
        result can be:
                        correct
                        undercall
                        overcall
                        miscalled_sib
                        miscall
        
        Stategy:
        First see if the NCBI tax IDs agree. If so, we're done, regardless of rank.
        If not, we'll need to use the taxonomy to figure it out. 
        """
        
        # Get our taxons
        true_tax = taxons.lookup_tax_id(ncbi_tax_id__true)
        true_taxonomy, true_taxonomy_ids, true_taxonomy_string = trim_and_unpack_taxon(true_tax, target_rank)
        
        test_tax = taxons.lookup_tax_id(ncbi_tax_id__test)
        test_taxonomy, test_taxonomy_ids, test_taxonomy_string = trim_and_unpack_taxon(test_tax, target_rank)

        # NCBI tax IDs?
        
        
        if test_taxonomy_ids[-1] == true_taxonomy_ids[-1]:
                return {
                        'result': 'correct',
                        'ranks_off': 0,
                        'target_rank': target_rank,
                        'lineage_true': true_taxonomy_string,
                        'lineage_classified': test_taxonomy_string,
                        'taxon_true': true_tax,
                        'taxon_classified': test_tax
                        }
                
        # If the IDs don't match, test the names. The true name is stored as organism. The test name is stored as name
        elif organism == name:
                return {
                        'result': 'correct',
                        'ranks_off': 0,
                        'target_rank': target_rank,
                        'lineage_true': true_taxonomy_string,
                        'lineage_classified': test_taxonomy_string,
                        'taxon_true': true_tax,
                        'taxon_classified': test_tax
                        }
        # Else. Use the taxonomy and figure it out from there
        
        
        # Check to see if we're a misscalled sib. Where the parent is the same
        elif len(test_taxonomy) > 1 and test_taxonomy_ids[-2] == true_taxonomy_ids[-2]:
                return {
                        'result': 'miscalled_sib',
                        'ranks_off': 1,
                        'target_rank': target_rank,
                        'lineage_true': true_taxonomy_string,
                        'lineage_classified': test_taxonomy_string,
                        'taxon_true': true_tax,
                        'taxon_classified': test_tax
                        }
                
        # Else, see if the end of our test taxonomy ID is in our true. Then were UNDERCALLED and we can quickly calculate by how much
        elif test_taxonomy_ids[-1] in true_taxonomy_ids:
                ranks_off = len(true_taxonomy_ids) - true_taxonomy_ids.index(test_taxonomy_ids[-1])-1
                return {
                        'result': 'undercall',
                        'ranks_off': ranks_off,
                        'target_rank': target_rank,
                        'lineage_true': true_taxonomy_string,
                        'lineage_classified': test_taxonomy_string,
                        'taxon_true': true_tax,
                        'taxon_classified': test_tax
                        }
                
        # Else. Now check if we're OVERCALLED, where our true tax ID is in our TEST lineage
        elif true_taxonomy_ids[-1] in test_taxonomy_ids:
                ranks_off = len(test_taxonomy_ids) - test_taxonomy_ids.index(true_taxonomy_ids[-1])-1
                return {
                        'result': 'overcall',
                        'ranks_off': ranks_off,
                        'target_rank': target_rank,
                        'lineage_true': true_taxonomy_string,
                        'lineage_classified': test_taxonomy_string,
                        'taxon_true': true_tax,
                        'taxon_classified': test_tax
                        }
        
        else:
        # Implicit Else. Must be a true miscall
                return {
                        'result': 'miscall',
                        'ranks_off': calculate_ranks_off(test_taxonomy_ids,true_taxonomy_ids,1),
                        'target_rank': target_rank,
                        'lineage_true': true_taxonomy_string,
                        'lineage_classified': test_taxonomy_string,
                        'taxon_true': true_tax,
                        'taxon_classified': test_tax
                        }
                
        
        

def test_row(row, taxons, target_rank=None):
        """
        For each row, try to figure out if the classification matches the true organism
        for a given target rank
        INPUTS:
                row: A row from a dataframe
                taxons: A TaxonCache instance
                target_rank: Rank at which to match. If none, try to match perfectly. Not implemented yet. 

        
        OUTPUT: A tuple of (test_result, ranks_off)
                test_result can be:
                        correct
                        undercalled
                        overcalled
                        miscalled_sib
                        miscall
                        dropped
        
        Issue to deal with:
        Some classifiers actually give a small LIST of possible taxIDs and names, essentially saying,
        "it's one of these". Ex Strep sangunins Or Strep parasanguinis.
        
        To deal with this, we will loop through EACH possible answer, and score them, weighting the answers.
        """
        num_classifications = len(row.ncbi_tax_id__test)
        # See if we need to loop at all
        if num_classifications == 1:
                return [test_match(taxons, row.ncbi_tax_id__true, row.ncbi_tax_id__test[0], row.organism, row.name__classified, target_rank)]
        else:
                classifications = []
                for ncbi_tax_id__test in row.ncbi_tax_id__test:
                        classifications.append(test_match(taxons, row.ncbi_tax_id__true, ncbi_tax_id__test, row.organism, row.name__classified, target_rank))
                
                return classifications
        
        
        
        
   

def main():
        args_parser = argparse.ArgumentParser()
        
        args_parser.add_argument('True', help='INPUT truth (map) file, in CSV format with the source of each sequence, by seq ID')
        args_parser.add_argument('Test', help='INPUT CSV file, listing by sequence ID the classifier\'s output')
        args_parser.add_argument('output_prefix', help='Prefix to use for the outputs')
        args_parser.add_argument('--strip','-s', help='Strip test seq IDs', action='store_true')
        args_parser.add_argument('--email','-e',required=True)
        args_parser.add_argument('--url','-u',help='Base URL to REST db with taxonomic information. If not provided, we will default to NCBI')
        #args_parser.add_argument('--rank','-R', help='Rank (species, genus, etc) to Target for SUMMARY by community. Default is none.', default=None, choices=(None,'species','genus'))
        
        
        
        args = args_parser.parse_args()
        
        
        if args.url:
            taxons = TaxonCache(source='custom', url = args.url, email=args.email)
        else:
            taxons = TaxonCache(email=args.email)
        
        # Load data
        true_df = pd.read_csv(args.True)
        test_df = pd.read_csv(args.Test)
        
        # We need to clean up the test_df a bit, unfortunately. Name is used as a colum, but reserved in pandas
        test_df['name__classified'] = test_df.name
        test_df.drop('name', axis=1, inplace=True)
        
        # IF our user has asked to strip the test seq IDs, do so here
        if args.strip:
                test_df.seq = test_df.seq.apply(lambda sid: sid.split('-')[0])
        
        # Similarly, we should convert our ncbi_tax_id to integers
        true_df.ncbi_tax_id = true_df.ncbi_tax_id.apply(lambda t: int(t))
        
        
        # Merge the frames
        #merged_df = pd.merge(true_df, test_df, how='inner', left_on='seqID', right_on='seq', suffixes=('__true','__test'))
        merged_df = pd.merge(true_df, test_df, how='outer', left_on='seqID', right_on='seq', suffixes=('__true','__test'))
        # Convert string lists of tax ids to tuples
        merged_df['ncbi_tax_id__test'] =merged_df.ncbi_tax_id__test.apply(handle_tax_ids)
        
        #merged_df['results_raw'] = merged_df.apply(lambda r: test_row(r,taxons), axis=1)
                
        
        json_data = {
                'taxons':               {},     # Dict of taxons. Key is NCBI tax ID. Value is a set of nested dicts for the taxon, including lineage
                'communities':          [],     # A list of community names
                'classifications':      [],     # A list of nested dicts. Each list item is:
                                                # 'community':          community name
                                                # 'source':             source sequence ID
                                                # 'num_calls_for_true': How many calls did the classifier make for this community and source?
                                                # 'name_true':          True name of source sequence (NOT restricted to target rank)
                                                # 'name_classified':    Classified name (NOT restricted to target rank)
                                                # 'ncbi_tax_id__true':  TRUE NCBI tax id (not limited to target rank). Lookup in taxons if you want...
                                                # 'ncbi_tax_id__test':  TEST NCBI tax id (not limited to target rank). Lookup in taxons if you want...
                                                # 'target_rank':        Rank targeted with this answer
                                                # 'lineage_true':       TRUE lineage, restricted to target rank, in form of semi-delimited string
                                                # 'lineage_test':       TEST lineage, restricted to target rank, in form of semi-delimited string
                                                # 'result':             RESULT of this classification at the target rank
                                                # 'ranks_off':          How many ranks off was this classification, given the target rank
                                                # 'weight':             Weight (in reads) for this call
        }
        
        json_fn = args.output_prefix+".json"
        json_f = open(json_fn,'w')
        
        community_f = {}
        community_w = {}
        for t_rank in TARGET_RANKS:
                if not t_rank:
                        t_rank = 'none'
                community_f[t_rank] = open(args.output_prefix+"."+t_rank+'.csv','w')
                community_w[t_rank] = csv.writer(community_f[t_rank])
                community_w[t_rank].writerow([
                                'community',
                                'correct',
                                'undercalled',
                                'overcalled',
                                'miscalled_sibs',
                                'miscalled',
                                'dropped',
                            ])
        
        
        for comm_name, comm_block in merged_df.groupby('community__true'):
        # Group by community, and start calculations here
            print "For Organisms in ", comm_name
            correct = {}
            undercall = {}
            overcall = {}
            misscall = {}
            misscalled_sib = {}
            dropped = {}
            
            for t_rank in TARGET_RANKS:
                if not t_rank:
                        t_rank = 'none'
                correct[t_rank]         = 0.0
                undercall[t_rank]       = 0.0
                overcall[t_rank]        = 0.0
                misscall[t_rank]        = 0.0
                misscalled_sib[t_rank]  = 0.0
                dropped[t_rank]         = 0.0
            
            json_data['communities'].append(comm_name)
            # From the perspective of TRUE distinct organisms
            
            for source_otu, source_block in comm_block.groupby('sourceSeq'):
                # For each true organism in this community, see how we did
                num_calls_for_true_otu = len(source_block.otu_id.unique())
                total_weight = float(len(source_block))
                
                
                # Test here for dropped true OTUs
                dropped_weight = len(source_block[pd.isnull(source_block.otu_id)]) / total_weight
                if dropped_weight > 0.0: # IF we have some mass that was dropped...
                        # Add the mass to our aggregate data
                        true_tax = taxons.lookup_tax_id(source_block.iloc[0].ncbi_tax_id__true)
                        json_data['taxons'][true_tax['TaxId']] = true_tax
                        lineage_true = trim_and_unpack_taxon(true_tax,None)[2]
                        # And then add an entry into our JSON table for this
                        for t_rank in TARGET_RANKS:
                                if not t_rank:
                                        t_rank = 'none'
                                dropped[t_rank] += dropped_weight
                                json_data['classifications'].append({
                                                                    'community':            comm_name,          #community name
                                                                    'source':               source_block.iloc[0].sourceSeq,  # source sequence ID
                                                                    'num_calls_for_true': num_calls_for_true_otu,  #How many calls did the classifier make for this community and source?
                                                                    'name_true':          source_block.iloc[0].organism,          # True name of source sequence (NOT restricted to target rank)
                                                                    'name_classified':    None, #Classified name (NOT restricted to target rank)
                                                                    'ncbi_tax_id__true':  true_tax['TaxId'],        # TRUE NCBI tax id (not limited to target rank). Lookup in taxons if you want...
                                                                    'ncbi_tax_id__test':  None,        #TEST NCBI tax id (not limited to target rank). Lookup in taxons if you want...
                                                                    'target_rank':        t_rank,               #Rank targeted with this answer
                                                                    'lineage_true':       lineage_true, # TRUE lineage, restricted to target rank, in form of semi-delimited string
                                                                    'lineage_test':       None, # TEST lineage, restricted to target rank, in form of semi-delimited string
                                                                    'result':             'dropped', # RESULT of this classification at the target rank
                                                                    'ranks_off':          None, # How many ranks off was this classification, given the target rank
                                                                    'weight':             dropped_weight,
                                                            })
                                      

                for otu_id, otu_block in source_block.groupby('otu_id'):
                    # Pick a representitive example from this OTU's block (they should all be the same)
                    rep_row = otu_block.iloc[0]
                    # Weight is a fraction of the number of reads with this call over the total number of reads for this source organism
                    call_weight = len(otu_block) / total_weight

                    # Look through the possible target ranks, only put out the summary for the one selected by our user:
                    for t_rank in TARGET_RANKS:
                                if not t_rank:
                                                t_rank = 'none'
                            
                                # Get our result.
                                testing_results = test_row(rep_row,taxons, t_rank) 
                                for testing_result in testing_results:
                                        result = testing_result['result']
                                            
                                            
                                        # Tick up our counter.
                                        if result == 'correct':
                                            correct[t_rank]+=call_weight
                                        elif result == 'undercall':
                                            undercall[t_rank]+=call_weight
                                        elif result == 'overcall':
                                            overcall[t_rank]+=call_weight
                                        elif result == 'miscalled_sib':
                                            misscalled_sib[t_rank]+=call_weight
                                        elif result == 'miscall':
                                            misscall[t_rank]+=call_weight
                                        else:
                                            print "Didn't understand result: ", result
                                            
                                        true_tax = testing_result['taxon_true']
                                        test_tax = testing_result['taxon_classified']
                                        
                                        # Update the taxons
                                        json_data['taxons'][true_tax['TaxId']] = true_tax
                                        json_data['taxons'][test_tax['TaxId']] = test_tax
                                        
                                        if not t_rank:
                                            t_rank = 'none'
                                        
                                        json_data['classifications'].append({
                                                        'community':            comm_name,          #community name
                                                        'source':               rep_row.sourceSeq,  # source sequence ID
                                                        'num_calls_for_true': num_calls_for_true_otu,  #How many calls did the classifier make for this community and source?
                                                        'name_true':          rep_row.organism,          # True name of source sequence (NOT restricted to target rank)
                                                        'name_classified':    rep_row.name__classified, #Classified name (NOT restricted to target rank)
                                                        'ncbi_tax_id__true':  true_tax['TaxId'],        # TRUE NCBI tax id (not limited to target rank). Lookup in taxons if you want...
                                                        'ncbi_tax_id__test':  test_tax['TaxId'],        #TEST NCBI tax id (not limited to target rank). Lookup in taxons if you want...
                                                        'target_rank':        t_rank,               #Rank targeted with this answer
                                                        'lineage_true':       testing_result['lineage_true'], # TRUE lineage, restricted to target rank, in form of semi-delimited string
                                                        'lineage_test':       testing_result['lineage_classified'], # TEST lineage, restricted to target rank, in form of semi-delimited string
                                                        'result':             result, # RESULT of this classification at the target rank
                                                        'ranks_off':          testing_result['ranks_off'], # How many ranks off was this classification, given the target rank
                                                        'weight':             call_weight,
                                                        
                                                })

                                
                        
        
            for t_rank in TARGET_RANKS:
                if not t_rank:
                        t_rank = 'none'
                print "Correct \t Undercalled \t Overcalled \t  Miscalled sibs \t Miscalled \t Dropped at "+t_rank
                print '\t', correct[t_rank],'\t|\t', undercall[t_rank], '\t|\t', overcall[t_rank],'\t|\t', misscalled_sib[t_rank], '\t|\t', misscall[t_rank], '\t|\t', dropped[t_rank]
                community_w[t_rank].writerow([
                    comm_name,
                    correct[t_rank],
                    undercall[t_rank],
                    overcall[t_rank],
                    misscalled_sib[t_rank],
                    misscall[t_rank],
                    dropped[t_rank]
                ])
            
            
        json.dump(json_data,json_f)
        json_f.close()
        for t_rank in community_f:
                community_f[t_rank].close()
        
                
                
        



if __name__ == '__main__':
    main()
