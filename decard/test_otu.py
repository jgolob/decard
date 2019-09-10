#!/usr/bin/env python
import argparse
import pandas as pd
import numpy as np
import csv

def main():
    args_parser = argparse.ArgumentParser()

    args_parser.add_argument('True')
    args_parser.add_argument('Test')
    args_parser.add_argument('output')

    args = args_parser.parse_args()
    
    true_df = pd.read_csv(args.True)
    test_df = pd.read_csv(args.Test)

    o_f = open(args.output,'w')
    writer = csv.writer(o_f)
    writer.writerow(['Community','Correct_Match','Correct_Split','Incorrect_Match','Incorrect_Split','Dropped','Total','sensitivity','specificity','ppv','npv','accuracy'])    
    merged_df = pd.merge(true_df, test_df, how='left', left_on='seqID', right_on='seq', suffixes=('__true','__test'))
    try:
        merged_df.set_index('seqID', inplace=True)
    except KeyError:
        merged_df.set_index('seqID__true', inplace=True)
        
    for comm_name, comm_block in merged_df.groupby('community__true'):
        print comm_name
        correct_match = 0
        correct_split = 0
        incorrect_match = 0
        incorrect_split = 0
        # Limit to not dropped
        comm_block_kept = comm_block[~comm_block.otu_id.isnull()]
        dropped = len(comm_block.index.unique()) - len(comm_block_kept.index.unique())
        for i in xrange(1,len(comm_block_kept)):
            left = comm_block_kept[i:]
            right = comm_block_kept[:-i]
            true_clusters = np.array(left.reset_index(drop=True).sourceSeq == right.reset_index(drop=True).sourceSeq)
            test_clusters = np.array(left.reset_index(drop=True).otu_id == right.reset_index(drop=True).otu_id)
            correct = (true_clusters == test_clusters)
            incorrect = ~correct
            correct_match += len([m for m in test_clusters[correct] if m])
            correct_split += len([m for m in test_clusters[correct] if not m])
            incorrect_match += len([m for m in test_clusters[incorrect] if m])
            incorrect_split += len([m for m in test_clusters[incorrect] if not m])
        print "Correct Match \t|\t Correct Split \t|\t Incorrect Match \t|\t Incorrect Split \t|\t Dropped \t|\t Total"
        print correct_match,'\t', correct_split,'\t' ,incorrect_match,'\t', incorrect_split, '\t', dropped, len(comm_block.index.unique())
        writer.writerow([comm_name,correct_match,correct_split,incorrect_match,incorrect_split,dropped,len(comm_block.index.unique()),
                         float(correct_match) / (correct_match+incorrect_split), # sensitivity
                         float(correct_split) / (correct_split+incorrect_match), # specificity
                         float(correct_match) / (correct_match+incorrect_match), # PPV
                         float(correct_split) / (correct_split+ incorrect_split), # NPV
                         (float(correct_match)+correct_split) / (correct_match+correct_split+incorrect_match+incorrect_split) # Accuracy
                         ])

    o_f.close()
    


if __name__ == '__main__':
    main()
