#!/usr/bin/env python
import argparse
import pandas as pd
import texttable
import sys
import os.path
import csv

def output_latex(result_table):
    import pylatex
    doc = pylatex.Document()
    """ TO BE IMPLEMENTED """
    
    
    
def normalize_results(result):
    result_n = pd.DataFrame()
    result_n['total'] = result.correct + result.undercalled + result.miscalled + result.miscalled_sibs + result.overcalled+result.dropped
    result_n['correct'] = result.correct / result_n.total*100.0
    result_n['undercalled'] = result.undercalled / result_n.total*100.0
    result_n['miscalled'] = result.miscalled / result_n.total*100.0
    result_n['miscalled_sibs'] = result.miscalled_sibs / result_n.total*100.0
    result_n['overcalled'] = result.overcalled / result_n.total*100.0
    result_n['dropped'] = result.dropped / result_n.total*100.0
    
    return result_n

def generate_result_table(classifications, names = None):
    resultTable = []
    

    
    if not names:
        names = ["Input "+str(i) for i in xrange(len(classifications))]
    
    for name, result in zip(names, classifications):
        resultTable.append({
            'name':     name,
            'correct':  {'mean': result.correct.mean(),
                         'std': result.correct.std(), },
            'undercalled':  {'mean': result.undercalled.mean(),
                             'std': result.undercalled.std(),},
            'miscalled_sibs':    {'mean': result.miscalled_sibs.mean(),
                             'std': result.miscalled_sibs.std(),},
            'overcalled':   {'mean': result.overcalled.mean(),
                             'std': result.overcalled.std(), },
            'miscalled':    {'mean':  result.miscalled.mean(),
                             'std':   result.miscalled.std(),},
            'dropped':    {'mean':  result.dropped.mean(),
                             'std':   result.dropped.std(),},
        })
    return resultTable

def generate_text_table(result_table):
    resultTable_textTable = texttable.Texttable(max_width=160)
    # Headers
    resultTable_textTable.header(['Input', 'Correct', "",'Undercalled',"", 'Miscalled Sibling',"", 'Overcalled',"", 'Miscalled', "", "Dropped", ""])
    resultTable_textTable.add_row(['', 'm', "std", 'm', "std", 'm', "std", 'm', "std", 'm', "std", 'm', "std"])
    
    resultTable_textTable.add_rows([
        [r['name'],
         r['correct']['mean'],r['correct']['std'],
         r['undercalled']['mean'],r['undercalled']['std'],
         r['miscalled_sibs']['mean'],r['miscalled_sibs']['std'],
         r['overcalled']['mean'],r['overcalled']['std'],
         r['miscalled']['mean'],r['miscalled']['std'],
         r['dropped']['mean'],r['dropped']['std'],
         ] for r in result_table], header=False)
    
    return resultTable_textTable


def output_results(result_table, file_h=sys.stdout):
    writer = csv.writer(file_h)
    
    writer.writerow(['Input', 'Correct', "",'Undercalled',"", 'Miscalled Sibling',"", 'Overcalled',"", 'Miscalled', "", "Dropped", ""])
    writer.writerow(['', 'm', "std", 'm', "std", 'm', "std", 'm', "std", 'm', "std", 'm', "std"])
    writer.writerows([
        [r['name'],
         r['correct']['mean'],r['correct']['std'],
         r['undercalled']['mean'],r['undercalled']['std'],
         r['miscalled_sibs']['mean'],r['miscalled_sibs']['std'],
         r['overcalled']['mean'],r['overcalled']['std'],
         r['miscalled']['mean'],r['miscalled']['std'],
         r['dropped']['mean'],r['dropped']['std'],
         ] for r in result_table])
    
    

def main():
    args_parser = argparse.ArgumentParser()
    
    args_parser.add_argument('inputs', nargs='+', help='Classification test files to compare')
    args_parser.add_argument('--normalize','-n', action='store_true', help='If selected, the data from inputs will be normalized to a percentage')
    args_parser.add_argument('--output','-o', help='Optional file into which to output the results. Else just to stdout')
    
    args = args_parser.parse_args()


    classifications = [pd.read_csv(filename) for filename in args.inputs]

    
    if args.normalize:
        classifications = [normalize_results(classification) for classification in classifications]
    
    names = [os.path.basename(fn) for fn in args.inputs]
    
    resultTable = generate_result_table(classifications,names)

    sys.stdout.writelines(generate_text_table(resultTable).draw())
    sys.stdout.write('\n')
    
    if args.output:
        with open(args.output,'w') as out_f:
            output_results(resultTable,out_f)
        
    
if __name__ == '__main__':
    main()   
        
    
    
