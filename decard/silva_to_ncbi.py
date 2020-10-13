#!/usr/bin/env python
import argparse
import sqlite3
import logging
import csv

logging.basicConfig(format='SILVA-NCBI:%(levelname)s:%(asctime)s: %(message)s', level=logging.INFO)

BANNED_NAMES = {
    'uncultured bacterium',
    'uncultured',
    'metagenome',
}

def main():
    args_parser = argparse.ArgumentParser(
        description="""Connect SILVA taxonomy to NCBI. 
        Needs a silva header and a taxonomy database made by taxtastic."""
    )

    args_parser.add_argument('--silva', '-s',
                             help="IN: Silva Headers",
                             required=True,
                             type=argparse.FileType('rt'))

    args_parser.add_argument('--taxdb', '-t',
                             help="IN: NCBI taxonomy database (SQLITE3) as made by taxtastic",
                             required=True)

    args_parser.add_argument('--ncbi', '-n',
                             help="OUT: Silva <-> NCBI in CSV format",
                             required=True,
                             type=argparse.FileType('wt'))

    args = args_parser.parse_args()
    
    logging.info("Loading taxonomy database")

    taxdb_con = sqlite3.connect(args.taxdb)
    taxdb_cur = taxdb_con.cursor()

    logging.info("Building database indicies")
    taxdb_cur.execute("CREATE INDEX IF NOT EXISTS names_taxid ON names (tax_id)")
    taxdb_cur.execute("CREATE INDEX IF NOT EXISTS names_name ON names (tax_name)")
    taxdb_cur.execute("CREATE INDEX IF NOT EXISTS nodes_taxid ON nodes (tax_id)")    

    logging.info("Loading SILVA Taxonomy")

    accession_NCBI = []
    ### Columns are 
    # silva_accession
    # silva lineage
    # ncbi_taxid
    # ncbi_taxname
    # ncbi_rank
    # to_species

    for entry in args.silva:
        ncbi_taxid = None
        accession = entry.split()[0]
        lineage_str = " ".join(entry.split()[1:])
        lineage_list = [t for t in lineage_str.split(';') if t not in BANNED_NAMES]
        if len(lineage_list[-1].split()) > 1:
            species = lineage_list.pop()
            taxdb_cur.execute("SELECT names.tax_id, names.tax_name, nodes.rank FROM names JOIN nodes ON names.tax_id=nodes.tax_id WHERE names.tax_name == ?", (species,))
            res = taxdb_cur.fetchone()
            if res is not None:
                ncbi_taxid = res[0]
                ncbi_taxname = res[1]
                ncbi_rank = res[2]
                accession_NCBI.append((
                    accession,
                    lineage_str,
                    ncbi_taxid,
                    ncbi_taxname,
                    ncbi_rank,
                    True
                ))
                continue
            genus = species.split()[0]
            taxdb_cur.execute("SELECT names.tax_id, names.tax_name, nodes.rank FROM names JOIN nodes ON names.tax_id=nodes.tax_id WHERE names.tax_name == ?", (genus,))
            res = taxdb_cur.fetchone()
            if res is not None:
                ncbi_taxid = res[0]
                ncbi_taxname = res[1]
                ncbi_rank = res[2]
                accession_NCBI.append((
                    accession,
                    lineage_str,
                    ncbi_taxid,
                    ncbi_taxname,
                    ncbi_rank,
                    False
                ))
                continue
        # Implicit else
        while len(lineage_list) > 0:
            t = lineage_list.pop()
            
            taxdb_cur.execute("SELECT names.tax_id, names.tax_name, nodes.rank FROM names JOIN nodes ON names.tax_id=nodes.tax_id WHERE names.tax_name == ?", (t,))
            res = taxdb_cur.fetchone()
            if res is not None:
                ncbi_taxid = res[0]
                ncbi_taxname = res[1]
                ncbi_rank = res[2]
                accession_NCBI.append((
                    accession,
                    lineage_str,
                    ncbi_taxid,
                    ncbi_taxname,
                    ncbi_rank,
                    False
                ))
                break
            else:
                taxdb_cur.execute("SELECT names.tax_id, names.tax_name, nodes.rank FROM names JOIN nodes ON names.tax_id=nodes.tax_id WHERE names.tax_name == ?", (t.lower(),))
                res = taxdb_cur.fetchone()
                if res is not None:
                    ncbi_taxid = res[0]
                    ncbi_taxname = res[1]
                    ncbi_rank = res[2]
                    accession_NCBI.append((
                        accession,
                        lineage_str,
                        ncbi_taxid,
                        ncbi_taxname,
                        ncbi_rank,
                        False
                    ))
                    break                
                
        if ncbi_taxid is None:
            logging.error(
                "COULD NOT FIND NCBI taxon for accession {}".format(accession)
            )

    logging.info("Completed parsing and matching of silva to NCBI")

    logging.info("Writing output")
    ncbi_w = csv.writer(args.ncbi)
    # Header
    ncbi_w.writerow([
        'silva_accession',
        'silva_lineage',
        'ncbi_taxid',
        'ncbi_taxname',
        'ncbi_rank',
        'to_species'
    ])
    ncbi_w.writerows(accession_NCBI)
    logging.info("DONE!")

if __name__ == '__main__':
    main()
