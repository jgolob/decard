#!/usr/bin/env python

import json
import pandas as pd
import os


def read_classification_json(fn):
    """
        Read in the classification data from a json file.
        Input: File name
        Output: the unpacked data
            data['taxons']: A dictionary with a taxon id as a key, and a taxon object as the value
            data['communities']: A LIST of community names
            data['classifications]: A list of classifications. Each member is a dict:
                [u'name_true',          # Actual organism name
                        u'weight',      # Weight of this call for this TRUE organims
                        u'target_rank', # Rank targeted when judging the classification
                        u'community',   # From which community was this classification
                        u'ncbi_tax_id__test',   # tax ID of the classified result
                        u'source',              # Source accession (equivalent to TRUE OTU id)
                        u'name_classified',     # Name the classifier returned
                        u'ncbi_tax_id__true',   # TRUE NCBI tax ID (to use in taxon lookup table)
                        u'result',              # How it scored. correct, undercall, overcall, miscall, miscalled_sib, dropped
                        u'ranks_off',           # How many ranks off
                        u'num_calls_for_true',  # How many different classifications were there for this organism from this community?
                        u'lineage_true',        # TRUE Lineage in semi-delimited format
                        u'lineage_test']        # TEST lineage in semi-delimited format
            
    """
    with open(fn) as f:
        classification_data = json.load(f)
        f.close()
        
    return classification_data


def summarize_classification(cd, target_rank ='species'):
    classifications = [c for c in cd['classifications'] if c['target_rank']==target_rank]
    classification_df = pd.DataFrame(classifications)
    
    total = 0.0
    miscall_n = 0.0
    miscall_r = 0.0
    resolution_n = 0.0
    resolution_r = 0.0
    dropped_n = 0.0
    correct_n = 0.0
    
    
    for source, true_block in classification_df.groupby('source'):
        for n, row in true_block.iterrows():
            total += row.weight
            if row.result == 'correct':
                correct_n+=row.weight
            elif row.result == 'dropped':
                dropped_n+=row.weight
            elif row.result =='undercall' or row.result =='overcall':
                resolution_n+=row.weight
                resolution_r+=row.weight*row.ranks_off
            elif row.result =='miscall' or row.result == 'miscalled_sib':
                miscall_n+=row.weight
                miscall_r+=row.weight*row.ranks_off
            else:
                print row.result
                
    correct_p = (correct_n / total)*100.0
    miscall_p = (miscall_n/total)*100.0
    resolution_p = (resolution_n / total) *100.0
    dropped_p = (dropped_n / total)*100.0
    
    print correct_p, miscall_p, resolution_p, dropped_p
    
def summarize_classification_by_clade(cd, target_rank ='species', clade_level=4):
    classifications = [c for c in cd['classifications'] if c['target_rank']==target_rank]
    classification_df = pd.DataFrame(classifications)
    
    # Get the clade grouper
    
    classification_df['clade'] = classification_df.lineage_true.apply(lambda ls: ls.split('; ')[clade_level])
    clade_results = []
    print "Clade", "Correct", "Miscalled", "Resolution", "Dropped", "Miscall rank off", "Resolution rank off"
    for clade, clade_block in classification_df.groupby('clade'):
    
        total = 0.0
        miscall_n = 0.0
        miscall_r = 0.0
        resolution_n = 0.0
        resolution_r = 0.0
        dropped_n = 0.0
        correct_n = 0.0
        
        
        for source, true_block in clade_block.groupby('source'):
            for n, row in true_block.iterrows():
                total += row.weight
                if row.result == 'correct':
                    correct_n+=row.weight
                elif row.result == 'dropped':
                    dropped_n+=row.weight
                elif row.result =='undercall' or row.result =='overcall':
                    resolution_n+=row.weight
                    resolution_r+=row.weight*row.ranks_off
                elif row.result =='miscall' or row.result == 'miscalled_sib':
                    miscall_n+=row.weight
                    miscall_r+=row.weight*row.ranks_off
                else:
                    print row.result
                    
        correct_p = (correct_n / total)*100.0
        miscall_p = (miscall_n/total)*100.0
        resolution_p = (resolution_n / total) *100.0
        dropped_p = (dropped_n / total)*100.0
        
        miscall_r_w = miscall_r / total
        resolution_r_w = resolution_r / total

              
        print clade, correct_p, miscall_p, resolution_p, dropped_p, miscall_r_w, resolution_r_w
        clade_results.append(
        {
            'clade':    clade,
            'correct_p':    correct_p,
            'miscall_p':    miscall_p,
            'resolution_p': resolution_p,
            'dropped_p':    dropped_p,
            'miscall_rw':   miscall_r_w,
            'resolution_rw':    resolution_r_w,
        })
    return clade_results

def main():
    path = '.'
    dir_fn = os.listdir(path)
    json_fn = [fn for fn in dir_fn if fn[-5:]=='.json']
    
    classifications = [read_classification_json(os.path.join(path, fn)) for fn in json_fn]
    
    results_byClade = [summarize_classification_by_clade(cd) for cd in classifications]
    
    # Now convert to longitudinal format for subsequent stats
    clades = set()
    for result in results_byClade:
        clades.update([r['clade'] for r in result])
    
    long_format = []
    for classifier, results in zip(json_fn, results_byClade):
        for result in results:
            result['classifier'] = classifier
            long_format.append(result)
        
    long_df = pd.DataFrame(long_format)
    
    # First a simple view
    means_across_classifiers = []
    for clade, clade_block in long_df.groupby('clade'):
        means_across_classifiers.append([
            clade,
            clade_block.correct_p.mean(),
            clade_block.miscall_p.mean(),
            clade_block.resolution_p.mean(),
            clade_block.dropped_p.mean(),
            clade_block.miscall_rw.mean(),
            clade_block.resolution_rw.mean(),
        ])

    means_df = pd.DataFrame(means_across_classifiers, columns = ['clade','correct','miscall','resolution','dropped', 'miscalled_ranks', 'resolution_ranks'])        
    