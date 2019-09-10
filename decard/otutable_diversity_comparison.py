#!/usr/bin/env python
"""
        Goal: Compare diversity between the real community (as defined by a target.csv) and as classified (output format as otu_table)
        Input: xxxx.otutable.csv
        
"""
import argparse
import pandas as pd
import rpy2
from rpy2.robjects import r as R
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri
pandas2ri.activate()
import re
from Bio import SeqIO
import numpy as np
import os

import matplotlib
matplotlib.use('Agg')
from scipy import stats
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.gridspec as gridspec

header_re = re.compile('description=\"(?P<description>[^\"]+)\" organism=\"(?P<organism>[^\"]+)\" taxonomy=\"(?P<taxonomy>[^\"]+)\" ncbi_tax_id=\"(?P<ncbi_tax_id>[^\"]+)\"')


def get_species(row):
    """
        Try to figure out if we have a species name
    """
    if len(row.taxonomy_list)> 6: # If our taxonomy is long enough, use the seventh position as our species name
        return row.taxonomy_list[6]
    elif row.organism != row.taxonomy_list[-1]: # Do we have an organism name different than our last spot in the taxonomy list? If so, use that
        return row.organism
    else:
        return "" # Nothing useful. 


def get_taxonomy_from_source(row, rRNA_dir):
        try:
                seq_recs = SeqIO.parse(os.path.join(rRNA_dir,row.source_file),'fasta')
        
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

def make_phyloseq_object(otu_matrix_df, taxon_df):
    # Load the requisite R libraries
    R_base = importr('base')
    R_phyloseq = importr('phyloseq')
    
    # Convert our pandas data frames to the appropriate R data types
    otu_matrix_R = pandas2ri.py2ri(otu_matrix_df)
    taxon_df_R = pandas2ri.py2ri(taxon_df)
    taxon_matrix_R = R_base.as_matrix(taxon_df_R)
    
    # Create our phyloseq object
    phyloseq_d = R_phyloseq.phyloseq(R_phyloseq.otu_table(otu_matrix_R, taxa_are_rows = True),R_phyloseq.tax_table(taxon_matrix_R))
    
    return phyloseq_d
            
def load_target(target_fn, rRNA_dir):
    # Load
    target_df = pd.read_csv(target_fn)
    
    # Convert weights into counts
    target_df['count'] = target_df.weight / target_df.weight.min()
    
    # Get our taxonomies
    target_df['taxonomy'] = target_df.apply(lambda r: get_taxonomy_from_source(r, rRNA_dir),axis=1)
    target_df['taxonomy_list'] = target_df.taxonomy.apply(lambda t: t.split(';'))
    
    # Start in on our matrix
        
    # Get our source IDS
    source_IDs = target_df.sequence_id.unique()
    
    # OTU matrix
    otu_matrix = pd.DataFrame(index=source_IDs, columns=target_df.community_name.unique())
    
    # Fill the matrix
    for source_id, source_block in target_df.groupby('sequence_id'):
            # Then get the total weight for each OTU, by community
            community_weights = source_block.groupby('community_name').sum()['count']
            for comm, weight in community_weights.iteritems():
                # Use the loc accessor to put in the value
                otu_matrix.loc[source_id,comm]= weight
    # Add zeros for the remainder
    otu_matrix.fillna(0.0, inplace=True)
    # Convert to integers
    otu_matrix = otu_matrix.applymap(np.int64)

    # Next the taxon list
    
    tax_table = pd.DataFrame(index=source_IDs, columns=['Domain', 'Phylum', 'Class', 'Order','Family','Genus', 'Species'])
    # This isn't very pandas, but it'll work for now
    for sequence_id, row in tax_table.iterrows():
            #print sequence_id, len(target_df[target_df.sequence_id==sequence_id])
            taxonomy_list = target_df[target_df.sequence_id==sequence_id].iloc[0].taxonomy_list
            species = target_df[target_df.sequence_id==sequence_id].iloc[0].species
            # Populate
            row.Domain = lookup_rank(taxonomy_list,0)
            row.Phylum = lookup_rank(taxonomy_list,1)
            row.Class = lookup_rank(taxonomy_list,2)
            row.Order = lookup_rank(taxonomy_list,3)
            row.Family = lookup_rank(taxonomy_list,4)
            row.Genus = lookup_rank(taxonomy_list,5)
            row.Species = species
            
    # Now assemble into a phyloseq object.
    return make_phyloseq_object(otu_matrix, tax_table)

def load_otutable(fn, seq_info_df):
    otu_table = pd.read_csv(fn)
    # Do a bit of normalizing here, first make sure that our OTU table has a weight column. If not, default to a weight of 1
    if not 'weight' in otu_table.columns:
        otu_table['weight'] = 1
        
    # Next be sure that we have a community column, and that it actually has community data in it.
    if not 'community' in otu_table.columns or len(otu_table.community[otu_table.community.isnull()]) > 0:
        # if we do not have seq info data, stop here.
        if not hasattr(seq_info_df,'to_json'):
            sys.stderr.write("No community information avaiable. Consider providing a seqinfo file")
            return -1
        else:
            # Drop the bad community column
            try:
                otu_table.drop('community', axis=1, inplace=True)
            except ValueError: # Wasn't one in the first place. Fine
                pass
            otu_table = pd.merge(otu_table,seq_info_df,on='seq')
    
    # Create the OTU matrix, with each community getting a column, each OTU getting a row.
    otu_matrix_df = pd.DataFrame(columns=otu_table.community.unique(), index=otu_table.otu_id.unique())

    # Fill our matrix by looping through, first grouping by OTU_ID
    for otu_id, otu_block in otu_table.groupby('otu_id'):
        # Then get the total weight for each OTU, by community
        community_weights = otu_block.groupby('community').sum().weight
        for comm, weight in community_weights.iteritems():
            # Use the loc accessor to put in the value
            otu_matrix_df.loc[otu_id,comm]= weight
            
    # Fill the remaining with zeros
    otu_matrix_df.fillna(0, inplace=True)
    # Convert to integers
    otu_matrix_df = otu_matrix_df.applymap(np.int64)
    
    # Now work on our tax table, where the index is the otu_id.
    
    # first get our uniqe otu_id / taxon strings
    otu_strings = otu_table.drop_duplicates(['otu_id','taxonomy_string','name'])
    
    # Next set up our taxon_df
    taxon_df = pd.DataFrame()
    taxon_df['otu_id'] = otu_strings.otu_id
    taxon_df['taxonomy_string'] = otu_strings.taxonomy_string
    taxon_df['organism'] = otu_strings['name']
    # Then make our OTU ID the index
    taxon_df.set_index('otu_id', inplace=True)
     
    # Split our strings into lists
    taxon_df['taxonomy_list'] = taxon_df.taxonomy_string.apply(lambda ts: ts.split('; '))
    
    taxon_df['Domain'] = taxon_df.taxonomy_list.apply(lambda tl: tl[0] if len(tl) > 0 else "")
    taxon_df['Phylum'] = taxon_df.taxonomy_list.apply(lambda tl: tl[1] if len(tl) > 1 else "")
    taxon_df['Class'] = taxon_df.taxonomy_list.apply(lambda tl: tl[2] if len(tl) > 2 else "")
    taxon_df['Order'] = taxon_df.taxonomy_list.apply(lambda tl: tl[3] if len(tl) > 3 else "")
    taxon_df['Family'] = taxon_df.taxonomy_list.apply(lambda tl: tl[4] if len(tl) > 4 else "")
    taxon_df['Genus'] = taxon_df.taxonomy_list.apply(lambda tl: tl[5] if len(tl) > 5 else "")
    taxon_df['Species'] = taxon_df.apply(get_species,axis=1)
    
    taxon_df.drop('taxonomy_string', inplace=True, axis=1)
    taxon_df.drop('taxonomy_list', inplace=True, axis=1)
    taxon_df.drop('organism', inplace=True, axis=1)
    
    return make_phyloseq_object(otu_matrix_df,taxon_df)

    
    

def get_richness(phyloseq_d):
    R_phyloseq = importr('phyloseq')
    richness  =R_phyloseq.estimate_richness(phyloseq_d, split = True)
    return pandas2ri.ri2py(richness)
    
    
def main():
        args_parser = argparse.ArgumentParser()
        
        args_parser.add_argument('Target', help='INPUT  (target) file, in CSV format')
        args_parser.add_argument('--rRNA16SSU','-r', help="Directory containing 16SSU rRNA", required=True)
        args_parser.add_argument('--otutables', help="OTU table(s) in csv format", nargs='+', required=True)
        args_parser.add_argument('--seqinfo', '-si', help='(optional) Headerless CSV file mapping sequence ID to community.')
        args_parser.add_argument('--figure', '-F', help='OUTPUT (optional) filename where to put boxplots in PDF format')
        args_parser.add_argument('--statsfile', '-S', help='OUTPUT (optional) filename where to stats in Excel format')
        
        args = args_parser.parse_args()

        # If we have it, load our seq info
        if args.seqinfo:
            seq_info_df = pd.read_csv(args.seqinfo, names=['seq','community'])
        else:
            seq_info_df = None
        
        otutable_fns = args.otutables
        
        
        target_ps = load_target(args.Target, args.rRNA16SSU)
        
        otutable_pss = [load_otutable(fn,seq_info_df) for fn in otutable_fns]
        
        
        # Diversity adventure
        target_diversity = get_richness(target_ps)
        target_shannon = target_diversity.Shannon
        
        
        otutable_diversity = [get_richness(ps) for ps in otutable_pss]
        otutable_shannon = [otd.Shannon for otd in otutable_diversity]
        
        shannon_combined_df = pd.DataFrame(index=target_shannon.index)
        shannon_combined_df['target'] = target_shannon
        shannon_normalized_df = pd.DataFrame(index=target_shannon.index)
        for name, o_t_s in zip(otutable_fns, otutable_shannon):
            shannon_combined_df[name]= o_t_s
            shannon_normalized_df[name] = o_t_s / target_shannon
        
        if args.figure:
            # Box-whiskers plots
            fig = plt.figure(0,figsize=(5,10))
            
            gs = gridspec.GridSpec(2,1,height_ratios=[3,1])
            bw_ax = plt.subplot(gs[0])
            plt.title("Ratio of Estimated Shannon Diversity to True")
            
            bw_ax.boxplot(np.array(shannon_normalized_df), vert = False)
            leg_ax = plt.subplot(gs[1])
            leg_ax.axis('off')
            
            fontsize=6
            legtext = ""
            for i, label in enumerate(otutable_fns):
                legtext += str(i+1)+'. '+label+'\n'
            leg_ax.text(0, 0, legtext, fontsize=fontsize)

            
            
            pdf = PdfPages(args.figure)
            pdf.savefig(fig)
            pdf.close()
        
        if args.statsfile:   
            # Do some statistics
            # First, by paired t-test does our estimate differ from true?
            pval_v_true = []
            for i, fn in enumerate(otutable_fns):
                res = stats.ttest_rel(shannon_combined_df.target, shannon_combined_df.ix[:,i+1])
                # Depending on the scipy version, the result comes back in different ways
                try:
                    pval = res.pvalue
                except AttributeError:
                    pval = res[1]
                print fn, pval
                pval_v_true.append({
                    'classifier': fn,
                    'pval':     pval,
                    'mean':     shannon_combined_df.ix[:,i+1].mean(),
                    'std':      shannon_combined_df.ix[:,i+1].std(),
                    'n':        len(shannon_combined_df.ix[:,i+1])
                })
            
            pval_v_true_df = pd.DataFrame(pval_v_true)
            
            # Then all of the pairwise pvalues between classifiers
            pval_mat = pd.DataFrame(index=otutable_fns, columns=otutable_fns)
            for i in xrange(0,len(otutable_fns)):
                for j in xrange(i+1,len(otutable_fns)):
                    res = stats.ttest_rel(shannon_normalized_df.ix[:,i], shannon_normalized_df.ix[:,j])
                    try:
                        pval = res.pvalue
                    except AttributeError:
                        pval = res[1]
                    print otutable_fns[i], otutable_fns[j], pval
                    pval_mat.iloc[j,i] = pval
                    pval_mat.iloc[i,j] = pval
                    
            # Output time
            xw = pd.ExcelWriter(args.statsfile)
            
            shannon_combined_df.to_excel(xw,sheet_name='Shannon')
            shannon_normalized_df.to_excel(xw,sheet_name='Shannon_Normalized')
            pval_v_true_df.to_excel(xw,sheet_name='pval_v_true')
            pval_mat.to_excel(xw,sheet_name='pval_mat')
            xw.save()
            
        
            
        
if __name__ == '__main__':
    main()
