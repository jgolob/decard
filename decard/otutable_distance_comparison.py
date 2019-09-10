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
from rpy2.robjects import pandas2ri, numpy2ri
pandas2ri.activate()
import re
from Bio import SeqIO, Phylo
import numpy as np
import os

import StringIO

from scipy import stats
import sys

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

def make_phyloseq_object(otu_matrix_df, taxon_df, tree_r):
    # Load the requisite R libraries
    R_base = importr('base')
    R_phyloseq = importr('phyloseq')
    
    # Convert our pandas data frames to the appropriate R data types
    otu_matrix_R = pandas2ri.py2ri(otu_matrix_df)
    taxon_df_R = pandas2ri.py2ri(taxon_df)
    taxon_matrix_R = R_base.as_matrix(taxon_df_R)
    
    # Create our phyloseq object
    phyloseq_d = R_phyloseq.phyloseq(R_phyloseq.otu_table(otu_matrix_R, taxa_are_rows = True),R_phyloseq.tax_table(taxon_matrix_R), tree_r)
    
    return phyloseq_d

def get_clade(clade, rank_name):
    if not rank_name in clade['children']:
        clade['children'][rank_name] = {
            'children': {},
            'members':  [],
        }
    return clade['children'][rank_name]

def tree_dicts_to_newick(clade):
    
    members = [str(m)+':1' for m in clade['members']]
    
    if  len(clade['children']) ==0:
        return ','.join(members)
    else:  # we have some children

        child_strings = [tree_dicts_to_newick(clade['children'][c]) for c in clade['children'] ]
    
        return '('+','.join(members+child_strings)+':1)'      

def tree_dicts_to_R_tree(tree_dicts):
    R_phytools = importr('phytools')
    R_ape = importr('ape')
    
    # Get into a newick-like string
    newick_s = tree_dicts_to_newick(tree_dicts)+';'
    
    # Make the tree
    tree_r = R_phytools.read_newick(text=newick_s)
    # collapse singles and return
    return R_ape.collapse_singles(tree_r)
    
    

def tax_table_to_tree(tax_table):
    root = {
        'children':  {},
        'members':  [],
    }

    for idx, row in tax_table.iterrows():
        # Start at the top
        clade = root
        # Domain
        if row.Domain:
            clade = get_clade(clade, row.Domain)
        # Phylum
        if row.Phylum:
            clade = get_clade(clade, row.Phylum)
        # Class
        if row.Class:
            clade = get_clade(clade, row.Class)
        # Order
        if row.Order:
            clade = get_clade(clade, row.Order)
        # Family
        if row.Family:
            clade = get_clade(clade, row.Family)
        # Genus
        if row.Genus:
            clade = get_clade(clade, row.Genus)
        # Species
        if row.Species:
            clade = get_clade(clade, row.Species)
            
        clade['members'].append(idx)
        
    return root

            
def load_target(target_fn, rRNA_dir, tree_fn = None):
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
    
    if not tree_fn:
        # now make a tree from our tax table
        dict_tree=tax_table_to_tree(tax_table)
        # Make this into an R tree
        tree_r =tree_dicts_to_R_tree(dict_tree)
    else: # try to load in the tree
        R_phytools = importr('phytools')
        R_ape = importr('ape')
        tree_r = R_phytools.read_newick(tree_fn)
        tree_r = R_ape.collapse_singles(tree_r)
        
            
    # Now assemble into a phyloseq object.
    return make_phyloseq_object(otu_matrix, tax_table, tree_r)

def load_otutable(fn, seq_info_df, tree_fn=None):
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
    
    # Finally convert otu_ids to strings for consistency with the trees we're loading.
    otu_table.otu_id = otu_table.otu_id.apply(lambda o: str(o))
    
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
    
                
    R_phytools = importr('phytools')
    R_ape = importr('ape')
    if not tree_fn:
        # now make a tree from our tax table
        sys.stderr.write("MAKING TAXONOMIC TREE FOR "+fn+"\n")
        dict_tree=tax_table_to_tree(taxon_df)
        # Make this into an R tree
        tree_r =tree_dicts_to_R_tree(dict_tree)
    else: # We have a tree fn. Load it
        # First see if our OTU names agree with the tree decorators. 
        otu_id_set = set(otu_table.otu_id.unique())
        tree = Phylo.read(tree_fn,'newick')
        nodes_names = [t.name for t in tree.find_elements() if t.name]
        
        if len(otu_id_set.difference(nodes_names)) > 0: # we have some otu_ids not in the tree
            # Maybe we have sequence IDs in the tree instead of OTU ids (ala mothur)
            if len(otu_table[otu_table.seq == nodes_names[0]]) ==1:
                for node in tree.find_elements():
                    if node.name:
                        otu_id = otu_table[otu_table.seq == node.name].otu_id.iloc[0]
                        if otu_id in otu_id_set: # Only relabel once
                            node.name = otu_id
                            otu_id_set.remove(otu_id)
                        
                        
                # now make tree_r
                newick_s = StringIO.StringIO()
                Phylo.write(tree,newick_s,'newick')
                tree_r = R_phytools.read_newick(text=newick_s.getvalue())
            else:
                # Probably won't work, but try it anyways
                tree_r = R_phytools.read_newick(tree_fn)
        
        else: # All is fine.
            tree_r = R_phytools.read_newick(tree_fn)

        
        tree_r = R_ape.collapse_singles(tree_r)
        
    # Now assemble into a phyloseq object.
    return make_phyloseq_object(otu_matrix_df, taxon_df, tree_r)

    
def get_wunifrac_distance(phyloseq_d):
    R_phyloseq = importr('phyloseq')
    R_base = importr('base')
        
    distances = R_phyloseq.UniFrac(phyloseq_d, weighted=True, normalized=True, fast=True, parallel=False)
    distance_mat = R_base.as_matrix(distances)
    distance_df = pd.DataFrame(numpy2ri.ri2py(distance_mat),
                               index=pandas2ri.ri2py(R_phyloseq.sample_names(phyloseq_d)),
                                columns=pandas2ri.ri2py(R_phyloseq.sample_names(phyloseq_d))
                            )
    return distance_df

def get_distance_methods(with_dpcoa=False):
    R_phyloseq = importr('phyloseq')
    raw_list = pandas2ri.ri2py(R_phyloseq.distanceMethodList)
    method_groups = [r for r in raw_list]
    method_names = []
    for g in method_groups:
        method_names = method_names+[n for n in g]
    # Remove the generic entry, if it exists. Silently move on if it doesn't
    try:
        method_names.remove('ANY')
    except ValueError:
        pass
    
    if not with_dpcoa:
        # Remove the dpcoa entry, if it exists. Silently move on if it doesn't
        try:
            method_names.remove('dpcoa')
        except ValueError:
            pass

    return method_names

def get_distance(phyloseq_d, dist_method):
    R_phyloseq = importr('phyloseq')
    R_base = importr('base')
        
    distances = R_phyloseq.distance(phyloseq_d, method=dist_method)
    distance_mat = R_base.as_matrix(distances)
    distance_df = pd.DataFrame(numpy2ri.ri2py(distance_mat),
                               index=pandas2ri.ri2py(R_phyloseq.sample_names(phyloseq_d)),
                                columns=pandas2ri.ri2py(R_phyloseq.sample_names(phyloseq_d))
                            )
    return distance_df
    
    
def main():
        args_parser = argparse.ArgumentParser()
        
        args_parser.add_argument('Target', help='INPUT  (target) file, in CSV format')
        args_parser.add_argument('--rRNA16SSU','-r', help="Directory containing 16SSU rRNA", required=True)
        args_parser.add_argument('target_tree', help="File with tree in newick format for the target sequences")
        args_parser.add_argument('--otutables', help="OTU table(s) in csv format", nargs='+')
        args_parser.add_argument('--config', '-c', help="CSV file with input files. Use this OR otutables")
        args_parser.add_argument('--seqinfo', '-si', help='(optional) Headerless CSV file mapping sequence ID to community.')
        args_parser.add_argument('--figure', '-F', help='OUTPUT (optional) filename where to put boxplots in PDF format')
        args_parser.add_argument('--distfile', '-D', help='OUTPUT (optional) filename where to put the distances in Excel format')
        args_parser.add_argument('--dpcoa', help="Do dpcoa distance", action='store_true')
        args_parser.add_argument('--taxonomy', help='Use a taxonomy instead of a phylogeny', action='store_true')
        args_parser.add_argument('--justthree', help="Just do the first three distance measures", action='store_true')
        
        args = args_parser.parse_args()

        if (args.otutables and args.config):
            sys.stderr.write("Please only provide a config file OR otutables listed on the command line")
        elif not (args.otutables or args.config):
            sys.stderr.write("Please only provide either a config file OR otutables listed on the command line")
        elif args.config:
            config_df = pd.read_csv(args.config)
            if args.taxonomy:
                config_df['tree'] = None
            
        else:
            # No config file. Assemble one
            config_df = pd.DataFrame(columns=['name', 'otutable','tree'])
            config_df.otutable = args.otutables
            config_df['name'] = args.otutables
            config_df['tree'] = None
            
        config_df['load_success'] = None
        config_df['phyloseq'] = None

        # If we have it, load our seq info
        if args.seqinfo:
            seq_info_df = pd.read_csv(args.seqinfo, names=['seq','community'])
        else:
            seq_info_df = None
        
        if args.figure:
            # Only import all this if we're doing figures
            import matplotlib
            matplotlib.use('Agg')
            import matplotlib.pyplot as plt
            from matplotlib.backends.backend_pdf import PdfPages
            import matplotlib.gridspec as gridspec
            matplotlib.rcParams.update({'font.size': 8})
            pdf = PdfPages(args.figure)
            figure=True
        else:
            figure=False

        
        if args.distfile:
            xw = pd.ExcelWriter(args.distfile)
            distfile = True
            
        else:
            distfile = False
        
        if args.taxonomy:
            target_ps = load_target(args.Target,args.rRNA16SSU,None)
        else:
            target_ps = load_target(args.Target,args.rRNA16SSU,args.target_tree)
        
        for idx, row in config_df.iterrows():
                print idx, row
            #try:
                config_df.loc[idx,'phyloseq'] = load_otutable(row.otutable, seq_info_df, row.tree)
                config_df.loc[idx,'load_success'] = True
            #except:
            #    sys.stderr.write(row.otutable+" failed to load\n")
            #    config_df.loc[idx,'load_success'] = False

        
        # Only keep those successfully loaded
        config_df = config_df[config_df.load_success==True]
                
        distance_methods = get_distance_methods(args.dpcoa)
        
        if args.justthree:
            distance_methods = distance_methods[:3]
        
        for dist_method in distance_methods:
            sys.stdout.write(dist_method+'\n')
            # Get the distances for this method
            
            target_dm = get_distance(target_ps,dist_method)
            
            otutable_dms = {}
            for idx, row in config_df.iterrows():
                try:
                    otutable_dms[row['name']] = get_distance(row.phyloseq, dist_method)
                except:
                    print "Failed ", dist_method
                    pass
            
            # OK, now sort our distance data frames by the names, so that they are the same for each
            target_dm.sort_index(inplace=True)
            target_dm.sort_index(inplace=True, axis=1)
                
            # now calculate ratio of ESTIMATED / TRUE pairwise distance
            # Flatten our target
            target_da = np.array(target_dm)
            # Limit to upper triangle values
            target_pd = target_da[np.triu_indices(len(target_dm),1)]
            t_indicies = np.triu_indices(len(target_dm),1)
            dist_ratios = pd.DataFrame()
            dist_combined = pd.DataFrame()
            dist_combined['com_1'] = [target_dm.index[idx_1] for idx_1 in t_indicies[0]]
            dist_combined['com_2'] = [target_dm.columns[idx_2] for idx_2 in t_indicies[1]]
            dist_combined['target'] = target_pd
            dist_ratios['com_1'] = [target_dm.index[idx_1] for idx_1 in t_indicies[0]]
            dist_ratios['com_2'] = [target_dm.columns[idx_2] for idx_2 in t_indicies[1]]
            dist_ratios['target'] = target_pd

            
            for idx, row in config_df.iterrows():
                o_dm = otutable_dms[row['name']]
                
                dist_combined[row['name']] = dist_combined.apply(lambda r: o_dm.loc[r.com_1, r.com_2], axis = 1)
                dist_ratios[row['name']]    = dist_combined.apply(lambda r: o_dm.loc[r.com_1, r.com_2] / r.target, axis = 1)

            # Get rid of the target column in dist_ratios
            
            dist_ratios.drop('target', axis=1, inplace=True)
            
            if figure:
                plt.clf()
                classifier_names = ['target']+list(config_df['name'])
                # Density plots
                fig = plt.figure(0,figsize=(4,2*len(classifier_names)))
                fig.suptitle(dist_method)
                main_ax = fig.add_subplot(111)
                main_ax.set_xlabel('True Distance')
                main_ax.set_ylabel('Estimated Distance')
                # Turn off axis lines and ticks of the big subplot
                main_ax.spines['top'].set_color('none')
                main_ax.spines['bottom'].set_color('none')
                main_ax.spines['left'].set_color('none')
                main_ax.spines['right'].set_color('none')
                main_ax.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')
                
                for i, cn in enumerate(classifier_names):
                    spearman_r, spearman_p = stats.spearmanr(dist_combined.target, dist_combined[cn])
                    pearson_r, pearson_p = stats.pearsonr(dist_combined.target, dist_combined[cn])
                    print cn, spearman_r*spearman_r, spearman_p, pearson_r*pearson_r, pearson_p
                    ax = fig.add_subplot(len(classifier_names),1,i+1)
                    
                    ax.set_title(cn, fontsize=8)
                    ax.hist2d(dist_combined.target, dist_combined[cn], bins=100, cmap=plt.cm.gray_r)
                    ax.text(np.percentile(dist_combined.target,2.5),np.percentile(dist_combined[cn],2.5),"Spearman R^2= "+str(np.round(spearman_r*spearman_r,decimals=3)), fontsize=6)



                pdf.savefig(fig)
            
            
            if distfile:        
                # Do some statistics
                # First, by paired t-test does our estimate differ from true?
                # And some bootstrapping to get 95% CI of spearman
                pval_v_true = []
                num_bootstraps = 1000
                for i, fn in enumerate(list(config_df['name'])):
                    # BOOTSTRAPPING SPEARMANS
                    spearmans = []
                    list_idx = np.arange(len(dist_combined.target))
                    method_data = dist_combined[fn]
                    for j in xrange(num_bootstraps):
                        bootstrap_list_idx = np.random.choice(list_idx,size=len(list_idx), replace=True)
                        spearman_r, spearman_p = stats.spearmanr([dist_combined.target[idx] for idx in bootstrap_list_idx], [method_data[idx] for idx in bootstrap_list_idx])
                        spearmans.append(spearman_r)
                    res = stats.ttest_rel(dist_combined.target, method_data)
                    # Depending on the scipy version, the result comes back in different ways
                    try:
                        pval = res.pvalue
                    except AttributeError:
                        pval = res[1]
                    print fn, pval
                    pval_v_true.append({
                        'classifier': fn,
                        'pval':     pval,
                        'mean':     method_data.mean(),
                        'std':      method_data.std(),
                        'n':        len(method_data),
                        'spearman_r_2.5':   np.percentile(spearmans,2.5),
                        'spearman_r':   np.percentile(spearmans,50),
                        'spearman_r_97.5':   np.percentile(spearmans,97.5),
                    })
                
                pval_v_true_df = pd.DataFrame(pval_v_true)
                
                # Then all of the pairwise pvalues between classifiers
                pval_mat = pd.DataFrame(index=config_df['name'], columns=config_df['name'])
                for i in xrange(0,len(config_df['name'])):
                    for j in xrange(i+1,len(config_df['name'])):
                        res = stats.ttest_rel(dist_ratios[config_df['name'].iloc[i]], dist_ratios[config_df['name'].iloc[j]]) 
                        try:
                            pval = res.pvalue
                        except AttributeError:
                            pval = res[1]
                        print config_df['name'].iloc[i], config_df['name'].iloc[j], pval
                        pval_mat.iloc[j,i] = pval
                        pval_mat.iloc[i,j] = pval
                        
                # Output time
                dist_combined.to_excel(xw,sheet_name=dist_method, index=False)
                dist_ratios.to_excel(xw,sheet_name=dist_method+"_ratio", index=False)
                pval_v_true_df.to_excel(xw,sheet_name=dist_method+'_pval_v_true')
                pval_mat.to_excel(xw,sheet_name=dist_method+'_pval_mat')
    
        
        if distfile:
            xw.save()
        
        if figure:
            pdf.close()
                
        
            
        
if __name__ == '__main__':
    main()
