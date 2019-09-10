#!/usr/bin/env python

import argparse
import csv


def main():
    args_parser = argparse.ArgumentParser()

    args_parser.add_argument('map_file',  help='Sequences map file that identifies each seq')
    args_parser.add_argument('seq_info',  help='(OUTPUT) Where to put the seq_info output file')
    args_parser.add_argument('--tsv','-t', action='store_true', help='Tab seperated text file rather than comma-seperated file output')
    
    args = args_parser.parse_args()
    
    
    with open(args.map_file) as map_f:
        map_r = csv.reader(map_f)
        header = map_r.next()
        with open(args.seq_info,'w') as seq_info_f:
            if args.tsv:
                seq_info_w = csv.writer(seq_info_f, delimiter='\t')
            else:
                seq_info_w = csv.writer(seq_info_f)
            for row in map_r:
                seq_info_w.writerow(row[:2])
            
            seq_info_f.close()
        map_f.close()
            
if __name__ == "__main__":
    main()
