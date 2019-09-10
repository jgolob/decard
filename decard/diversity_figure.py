#!/usr/bin/env python

import argparse
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


def main():
    args_parser = argparse.ArgumentParser()
    args_parser.add_argument('input', help='Excel File containing diversity data')
    args_parser.add_argument('output_pdf', help='Filename where to put output (PDF format)')
    
    args = args_parser.parse_args()
    
    shannon = pd.read_excel(args.input, sheet='Shannon')
    classifier_names = shannon.columns
    
    fig = plt.figure(0,figsize=(4,2*len(classifier_names)))
    main_ax = fig.add_subplot(111)
    main_ax.set_xlabel('True Diversity')
    main_ax.set_ylabel('Estimated Diversity')
    # Turn off axis lines and ticks of the big subplot
    main_ax.spines['top'].set_color('none')
    main_ax.spines['bottom'].set_color('none')
    main_ax.spines['left'].set_color('none')
    main_ax.spines['right'].set_color('none')
    main_ax.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')
    
    
    for i, cn in enumerate(classifier_names):
        spearman_r, spearman_p = stats.spearmanr(shannon.target, shannon[cn])
        pearson_r, pearson_p = stats.pearsonr(shannon.target, shannon[cn])
        print cn, spearman_r*spearman_r, spearman_p, pearson_r*pearson_r, pearson_p
        ax = fig.add_subplot(len(classifier_names),1,i+1)
        ax.set_title(cn, fontsize=8)
        ax.scatter(shannon.target, shannon[cn])
        ax.text(shannon.target.min(),shannon[cn].max(),"Spearman R^2= "+str(spearman_r*spearman_r), fontsize=6)
        
    fig.tight_layout()
    
    pdf = PdfPages(args.output_pdf)
    pdf.savefig(fig)
    pdf.close()
    


if __name__ == '__main__':
    main()
