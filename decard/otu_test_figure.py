#!/usr/bin/env python
import argparse
import matplotlib

#matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
from matplotlib.backends.backend_pdf import PdfPages
from scipy.stats import ttest_rel



def otu_box_whiskers(source_dir, out_fn):
    otu_test_fns = os.listdir(source_dir)
    otu_test_fns.sort()
    otu_test_fns.reverse()
    otu_test_fns = [fn for fn in otu_test_fns if fn[-4:]=='.csv']
    
    otu_tests = []
    for otu_test_fn in otu_test_fns:
        otu_test = pd.read_csv(os.path.join(source_dir,otu_test_fn))
        otu_test['raw_name'] = otu_test_fn
        otu_test['percent_dropped'] = otu_test.apply(lambda r: float(r.Dropped) / float(r.Total)*100.0, axis = 1)
        otu_test['matthews'] = otu_test.apply(lambda r: (r.Correct_Match*r.Correct_Split - r.Incorrect_Match*r.Incorrect_Split)/np.sqrt(float(r.Correct_Match+r.Incorrect_Match)*(r.Correct_Match+r.Incorrect_Split)*(r.Correct_Split+r.Incorrect_Match)*(r.Correct_Split+r.Incorrect_Split)), axis=1)
        
        otu_tests.append(otu_test)
        
        

    
    specificities = [ot.specificity for ot in otu_tests]
    sensitivities = [ot.sensitivity for ot in otu_tests]
    matthews = [ot.matthews for ot in otu_tests]
    dropped = [ot.percent_dropped for ot in otu_tests]
    labels = [ot.iloc[0].raw_name[:-12] for ot in otu_tests]
    
    sensitivity_pvals = []
    specificity_pvals = []
    matthews_pvals = []
    dropped_pvals = []
    table_cols = []
    table_rows = []
    for i in xrange(0,len(otu_tests)):
        """
        sensitivity_row = [""]*(i+1)
        specificity_row = [""]*(i+1)
        matthews_row = [""]*(i+1)
        dropped_row = [""]*(i+1)
        """
        sensitivity_row = []
        specificity_row = []
        matthews_row = []
        dropped_row = []
        

        for j in xrange(0,len(otu_tests)):
            #print i,j, otu_tests[i].raw_name.iloc[0].split("_")[0], otu_tests[j].raw_name.iloc[0].split("_")[0]
            # Merge on community to ensure our pairs are correct
            if i!= j:
                
                merged = pd.merge(otu_tests[i], otu_tests[j], on='Community', suffixes=['__i', '__j'])
                sensitivity_pval = ttest_rel(merged.sensitivity__i, merged.sensitivity__j).pvalue
                specificity_pval = ttest_rel(merged.specificity__i, merged.specificity__j).pvalue
                matthews_pval = ttest_rel(merged.matthews__i, merged.matthews__j).pvalue
                dropped_pval  = ttest_rel(merged.Dropped__i, merged.Dropped__j).pvalue
                
                if sensitivity_pval < 0.01:
                    sensitivity_row.append('< 0.01')
                else:
                    sensitivity_row.append('%1.2f' % sensitivity_pval)
                
                if specificity_pval < 0.01:
                    specificity_row.append('< 0.01')
                else:
                    specificity_row.append('%1.2f' % specificity_pval)
                
                if matthews_pval < 0.01:
                    matthews_row.append('< 0.01')
                else:
                    matthews_row.append('%1.2f' % matthews_pval)
                
                if dropped_pval < 0.01:
                    dropped_row.append('< 0.01')
                else:
                    dropped_row.append('%1.2f' % dropped_pval)
                
            else:
                sensitivity_row.append('')
                specificity_row.append('')
                matthews_row.append("")
                dropped_row.append("")
            if i == 0:
                table_rows.append(otu_tests[j].iloc[0].raw_name[:-12])

        sensitivity_pvals.append(sensitivity_row)
        specificity_pvals.append(specificity_row)
        matthews_pvals.append(matthews_row)
        dropped_pvals.append(dropped_row)
        table_cols.append(otu_tests[i].iloc[0].raw_name.split("_")[0])
    
    
    # --
    fig= plt.figure(0, figsize=(11,8))
    spec_ax = plt.subplot(241)
    spec_ax.boxplot(specificities, labels=labels, vert=False)
    spec_ax.set_title('Specificity')
    xlocs, xlabels = plt.xticks()
    plt.setp(xlabels, rotation=90, fontsize=7)
    ylocs,ylabels = plt.yticks()
    plt.setp(ylabels, fontsize=7)
    
    sens_ax = plt.subplot(242)
    sens_ax.boxplot(sensitivities, labels=['']*len(labels),vert=False)
    xlocs, xlabels = plt.xticks()
    plt.setp(xlabels, rotation=90, fontsize=7)
    sens_ax.set_title('Sensitivity')
    
    matthews_ax = plt.subplot(243)
    matthews_ax.boxplot(matthews, labels=['']*len(labels),vert=False)
    xlocs, xlabels = plt.xticks()
    plt.setp(xlabels, rotation=90, fontsize=7)
    matthews_ax.set_title('Matthews')
    
    dropped_ax = plt.subplot(244)
    dropped_ax.boxplot(dropped, labels=['']*len(labels),vert=False)
    xlocs, xlabels = plt.xticks()
    plt.setp(xlabels, rotation=90, fontsize=7)
    dropped_ax.set_title('Dropped')
    
    spec_pval_table = plt.subplot(245)
    spec_pval_table.axis('tight')
    spec_pval_table.axis('off')
    spec_pval_table.table(
        cellText = specificity_pvals,
        rowLabels = table_rows,
        colLabels = table_cols,
        rowLoc='right',
        loc='center'
    )
    
    sens_pval_table = plt.subplot(246)
    sens_pval_table.axis('tight')
    sens_pval_table.axis('off')
    sens_pval_table.table(
        cellText = sensitivity_pvals,
        #rowLabels = labels[1:],
        colLabels = table_cols,
        loc='center'
    )
    
    matthews_pval_table = plt.subplot(247)
    matthews_pval_table.axis('tight')
    matthews_pval_table.axis('off')
    matthews_pval_table.table(
        cellText = matthews_pvals,
        #rowLabels = labels[1:],
        colLabels = table_cols,
        loc='center'
    )
    
    dropped_pval_table = plt.subplot(248)
    dropped_pval_table.axis('tight')
    dropped_pval_table.axis('off')
    dropped_pval_table.table(
        cellText = dropped_pvals,
        #rowLabels = labels[1:],
        colLabels = table_cols,
        loc='center'
    )
    
    fig.subplots_adjust(left=0.30)
    
    pdf = PdfPages(out_fn)
    pdf.savefig(fig)
    pdf.close()


def main():
    
    args_parser = argparse.ArgumentParser()
    
    args_parser.add_argument('otu_test_dir', help='Directory containing OTU test CSV files')
    args_parser.add_argument('output_pdf', help='Filename where to put output (PDF format)')
    
    args = args_parser.parse_args()
    
    otu_box_whiskers(args.otu_test_dir,args.output_pdf)
    
    

if __name__ == '__main__':
    main()