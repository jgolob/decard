#!/usr/bin/env python
import argparse
import os
import re
import subprocess
import time

def main():
    args_parser = argparse.ArgumentParser()


    args_parser.add_argument('otutable_dir', help="Directory with OTU table files")
    args_parser.add_argument('output_dir', help='Directory into which to put the outputs')
    args_parser.add_argument('True', help='INPUT truth (map) file, in CSV format with the source of each sequence, by seq ID')
    args_parser.add_argument('--strip','-s', help='Strip test seq IDs', action='store_true')
    args_parser.add_argument('--email','-e',required=True)
    args_parser.add_argument('--url','-u',help='Base URL to REST db with taxonomic information. If not provided, we will default to NCBI')
    args_parser.add_argument('--taskmanager','-t', help='Use DRMAA taskmanager to run tests in parallel', action='store_true')
    args_parser.add_argument('--no_otu', help='Skip OTU testing', action='store_true')
    args_parser.add_argument('--no_classification', help='Skip classification testing', action='store_true')
            
    args = args_parser.parse_args()
    
    if (args.no_otu and args.no_classification):
        print "No testing to be done. Exiting"
        return 0
    
    if args.taskmanager:
        from Drmaa_TaskManager import Drmaa_TaskManager
        taskmanager = Drmaa_TaskManager()
    
    # Check / make the output directories
    try:
        os.stat(args.output_dir)
    except:
        os.mkdir(args,output_dir)
    
    if not args.no_classification:
        try:
            os.stat(os.path.join(args.output_dir,'classification'))
        except:
            os.mkdir(os.path.join(args.output_dir,'classification'))
    
    if not args.no_otu:    
        try:
            os.stat(os.path.join(args.output_dir,'otu'))
        except:
            os.mkdir(os.path.join(args.output_dir,'otu'))
    
    re_otutable = re.compile('(?P<prefix>.+).otutable.csv')
    
    source_dir_files = os.listdir(args.otutable_dir)
    
    otutable_fn = [f for f in source_dir_files if re_otutable.match(f)]
    
    for otutable_f in otutable_fn:
        print otutable_f
        if not args.taskmanager:
            if not args.no_classification:
                if args.strip:
                    subprocess.call(['test_classification.py',
                                     args.True,
                                     os.path.join(args.otutable_dir,otutable_f),
                                     os.path.join(args.output_dir,'classification',re_otutable.match(otutable_f).group('prefix')),
                                     '-s',
                                     '-e '+args.email,
                                     '-u '+args.url,
                                    ])
                else:
                    subprocess.call(['test_classification.py',
                                     args.True,
                                     os.path.join(args.otutable_dir,otutable_f),
                                     os.path.join(args.output_dir,'classification',re_otutable.match(otutable_f).group('prefix')),
                                     '-e '+args.email,
                                     '-u '+args.url,
                                    ])
            if not args.no_otu:
                subprocess.call([
                    'test_otu.py',
                    args.True,
                    os.path.join(args.otutable_dir, otutable_f),
                    os.path.join(args.output_dir,'otu',re_otutable.match(otutable_f).group('prefix')+'.otutest.csv'),            
                    ])
        else: # Use taskmanager
            if not args.no_classification:
                if args.strip:
                    taskmanager.startJob('test_classification.py',
                                     parameters = [args.True,
                                     os.path.join(args.otutable_dir,otutable_f),
                                     os.path.join(args.output_dir,'classification',re_otutable.match(otutable_f).group('prefix')),
                                     '-s',
                                     '-e '+args.email,
                                     '-u '+args.url,
                                    ],
                                    num_cpus = 1,
                                    partition = 'campus,boneyard',
                                    time='0:15'
                                    )
                else:
                    taskmanager.startJob('test_classification.py',
                                     parameters = [args.True,
                                     os.path.join(args.otutable_dir,otutable_f),
                                     os.path.join(args.output_dir,'classification',re_otutable.match(otutable_f).group('prefix')),
                                     '-e '+args.email,
                                     '-u '+args.url,
                                    ],
                                    num_cpus = 1,
                                    partition = 'campus,boneyard',
                                    time='0:15'
                                    )
            
            if not args.no_otu:
                taskmanager.startJob('test_otu.py',
                                     parameters = [args.True,
                                    os.path.join(args.otutable_dir, otutable_f),
                                os.path.join(args.output_dir,'otu',re_otutable.match(otutable_f).group('prefix')+'.otutest.csv'),            
                                ],
                                    num_cpus = 1,
                                    partition = 'campus,boneyard',
                                    time='0:30'
                                    )
    
    if args.taskmanager:
        while taskmanager.jobs_current:
            print(chr(27) + "[2J")
            print "Status: "
            taskmanager.printStatus()
            time.sleep(10)
            

if __name__ == '__main__':
    main()
