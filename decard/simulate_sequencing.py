#!/usr/bin/env python
'''
simulate_sequencing.py
Derived from fasta2fastq.py by Tom Smith from the CGAT package, with respect to the initial authors.

=====================================

:Author: Jonathan Golob

Purpose
-------

This module takes a fasta file of simulated amplicons, and simulates sequencing of these reads.
This involves
-> Adding sequencing errors / uncertainty to the reads
-> (optionally) Converting to paired-end reads (for Illuminia-style reads)
-> (optionally) Convert to FASTQ files. 


Options
-------

--output-paired-end
   generate paired-end reads (defaults to single end)

--sequence-error-phred
   the sequencing error rate (phred scale)

--output-read-length
   the length of the outputted reads

--output-counts
   filename for counts per fasta entry

--output-quality-format
   the format of the sequence qualities (+33 = Sanger)

--insert-length-mean
   the mean insert length

--insert-length-sd
   the standard deviation for the insert length

--premrna-fraction
   the fraction of reads to simulate from pre-mRNA. Default is 0.
   If set, must provide a pre-mRNA fasta file with:
      --infile-premrna-fasta

If generating paired end reads, the second outfile must be specified with:

--output-fastq2


Usage
-----
Type::

   python fasta2fastq.py --help

for command line help.



'''
import sys
import random
import numpy as np
import collections
import argparse
from Bio import SeqIO


def addSeqErrors(read=None, error_rate=10):
    ''' add sequencing errors to a read.
    Error rates are Phred scaled, so 30 = 1/1000'''

    error_rate = 10**(error_rate/-10.0)

    errors_dict = {"G": ["C", "T", "A"],
                   "C": ["G", "T", "A"],
                   "T": ["C", "G", "A"],
                   "A": ["C", "T", "G"],
                   "N": ["C", "T", "G", "A"]}

    probs = np.random.rand(len(read))
    return "".join([base if prob > error_rate and base != "N"
                    else random.choice(errors_dict[base])
                    for prob, base in zip(probs, read)])


def reverseComp(seq):
    ''' return the reverse complement sequence '''

    comp = {"G": "C",
            "C": "G",
            "A": "T",
            "T": "A",
            "N": "N"}

    return "".join([comp[base] for base in seq[::-1]])


def generateRead(entry, read_length=50, error_rate=40, paired=False,
                 insert_mean=0, insert_sd=1):
    ''' generate a read (or read pair) at random from a fasta entry for
    the given read length with sequencing errors according to error
    rate'''

    if paired:

        position = "not_OK"

        while position != "OK":

            r1_start = random.randint(0, len(entry)-read_length)
            r2_start = (r1_start + read_length +
                        int(np.random.normal(insert_mean, insert_sd)))

            if (r2_start <= (len(entry) - read_length) and r2_start >= r1_start):

                position = "OK"

                read1 = entry[r1_start: r1_start+read_length]
                read2 = reverseComp(
                    entry[r2_start: r2_start+read_length])

                final_read1 = addSeqErrors(read1, error_rate)
                final_read2 = addSeqErrors(read2, error_rate)

                return final_read1, final_read2

    else:
        start = random.randint(0, len(entry)-read_length)
        read = entry[start:start+read_length]

        final_read = addSeqErrors(read, error_rate)

        return final_read


def main():
    parser = argparse.ArgumentParser()

    
    parser.add_argument(
        "--output-quality-format", dest="q_format", type="int",
        help="sequence quality format, e.g 33 = +33/Sanger"
        "[default=%default].")

    parser.add_option(
        "--output-paired-end", dest="paired", action="store_true",
        help="generate paired end reads [default = %default].")

    parser.add_option(
        "--insert-length-mean", dest="insert_mean", type="float",
        help="mean insert length [default = %default].")

    parser.add_option(
        "--insert-length-sd", dest="insert_sd", type="float",
        help="insert length standard deviation [default = %default].")

    
    parser.add_option(
        "--output-read-length", dest="read_length", type="int",
        help="read length [default = %default].")

    parser.add_option(
        "--sequence-error-phred", dest="phred", type="int",
        help="phred quality score [default = %default].")

    parser.add_option(
        "--output-counts", dest="output_counts", type="string",
        help="name for counts outfile [default=%default].")

    parser.add_option(
        "--output-fastq2", dest="fastq2_out", type="string",
        help="filename for second fastq outfile [default=%default].")

    parser.add_option(
        "--infile-fasta", dest="infile_fasta", type="string",
        help="filename for simulated amplicon fasta[default=%default].")

    parser.set_defaults(
        q_format=33,
        paired=False,
        insert_mean=0,
        insert_sd=1,
        read_length=50,
        fastq2_out=None,
        output_counts=None,
        phred=30,
        infile_fasta=None
    )

    options = parser.parse_args()

    if options.paired:
        assert options.fastq2_out, ("must specify a second fastq outfile for "
                                    "paired end (--output-fastq2)")
        outf2 = open(options.fastq2_out, "w")

    
    assert options.infile_fasta, ("must specfify the location of the"
                                       "fasta file for the amplicons")

    # the sequence quality string will always be the same so define here
    sequence_quality = chr(options.q_format + options.phred)
    qual = "".join([sequence_quality] * options.read_length)

    iterator = SeqIO.parse(options.infile_fasta,'fasta')

    # set a cut off of twice the read/pair length for short entries
    if options.paired:
        minimum_entry_length = (
            2 * ((options.read_length * 2) + options.insert_mean))
    else:
        minimum_entry_length = 2 * options.read_length

    c = collections.Counter()
    counts = collections.Counter()
    copies = collections.Counter()

    for seq_rec in iterator:


        # reject short fasta entries
        if len(seq_rec) < minimum_entry_length:
            E.info("skipping short amplicon: %s length=%i"
                   % (entry.title, len(entry.sequence)))
            c['skipped'] += 1
            continue

        else:
            c['not_skipped'] += 1

        if options.paired:
            fragment_length = (
                (2 * options.read_length) + options.insert_mean)
        else:
            fragment_length = options.read_length

        if "N" in seq_rec.seq.upper():
            print "fasta entry %s contains unknown bases ('N')" % entry_id

        read = generateRead(entry=entry.sequence.upper(),
                                read_length=options.read_length,
                                error_rate=options.phred,
                                paired=options.paired,
                                insert_mean=options.insert_mean,
                                insert_sd=options.insert_sd)

            if options.paired:
                r1, r2 = read
                h1 = "@%s_%i/1" % (entry_id, i)
                h2 = "@%s_%i/2" % (entry_id, i)

                options.stdout.write("\n".join((h1, r1, "+", qual)) + "\n")
                outf2.write("\n".join((h2, r2, "+", qual)) + "\n")

            else:
                h = "@%s_%i/1" % (entry_id, i)

                options.stdout.write("\n".join((h, read, "+", qual)) + "\n")

        if options.premrna_fraction:
            c['pre_counts'] += n_reads_pre
            c['pre_copies'] += n_copies_pre

            for i in range(0, n_reads_pre):

                read = generateRead(entry=pre_entry.sequence.upper(),
                                    read_length=options.read_length,
                                    error_rate=options.phred,
                                    paired=options.paired,
                                    insert_mean=options.insert_mean,
                                    insert_sd=options.insert_sd)

                if options.paired:
                    r1, r2 = read
                    h1 = "@%s_pre-mRNA_%i/1" % (entry_id, i)
                    h2 = "@%s_pre-mRNA_%i/2" % (entry_id, i)

                    options.stdout.write("\n".join((h1, r1, "+", qual)) + "\n")
                    outf2.write("\n".join((h2, r2, "+", qual)) + "\n")

                else:
                    h = "@%s_pre-mRNA_%i/1" % (entry_id, i)

                    options.stdout.write("\n".join((h, read, "+", qual)) + "\n")

    if options.paired:
        outf2.close()

    with IOTools.openFile(options.output_counts, "w") as counts_out:

        counts_out.write("%s\n" % "\t".join(("id", "read_count", "tpm")))

        sum_copies = sum(copies.values())
        sum_counts = sum(counts.values())

        for entry_id, count in counts.iteritems():
            tpm = 1000000 * (float(copies[entry_id]) / sum_copies)
            counts_out.write(
                "%s\n" % "\t".join(map(str, (entry_id, count, tpm))))

    E.info("Reads simulated for %i fasta entries, %i entries skipped"
           % (c['not_skipped'], c['skipped']))

    E.info("Simulated: %i reads (%i mRNA, %i pre-mRNA), "
           "%f transcripts (%f mRNA, %f pre-mRNA)" % (
               sum_counts + c['pre_counts'], sum_counts, c['pre_counts'],
               sum_copies + c['pre_copies'], sum_copies, c['pre_copies']))

    E.Stop()

if __name__ == "__main__":
    main()

