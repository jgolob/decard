#!/usr/bin/env python

import argparse
import csv


def main():
    args_parser = argparse.ArgumentParser()

    args_parser.add_argument('map_file',  help='Sequences map file that identifies each seq')
    args_parser.add_argument('labels',  help='(OUTPUT) Where to put the labels file')
    args_parser.add_argument('--suffix','-s', help='(optional) suffix to add to labels')
    
    args = args_parser.parse_args()
    if args.suffix:
        suffix = str(args.suffix)
    else:
        suffix=""
    
    completed = set()
    with open(args.map_file) as map_f:
        map_r = csv.reader(map_f)
        header = map_r.next()
        with open(args.labels,'w') as labels_f:
            labels_w = csv.writer(labels_f)
            labels_w.writerow(['specimen','labels'])
            for row in map_r:
                if not row[1] in completed:
                    labels_w.writerow([row[1],row[1]+suffix])
                    completed.add(row[1])
            
            labels_f.close()
        map_f.close()
            
if __name__ == "__main__":
    main()
