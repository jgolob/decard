# README #

Welcome to DECARD (Detailed Evaluation Creation and Analysis of Read Data).

### What is DECARD? ###
DECARD can be used to simulate amplicon-based microbiome experiments, and test how classification software performs.
DECARD can help with each step of the process:

1. Downloading references sequences from NCBI (NCBI_16SMicrobial_update.py) and SILVA (SILVA_selectFromGenera.py)
2. Curating the reference sequences by GENUS and SPECIES. (selectFromGenera.py)
3. Performing quality control on the reference (genera_qc.py)
4. Generate targets for a set of communities (generate_targets.py)
5. Perform in-silico PCR and generate amplicon sequences (generate_sequences.py)
6. (Optional). Perform in-silico sequencing of the amplicons, with error addition using something like ART

These amplicons or sequences can then be classified. The resultant classifications can be evaluated
1. For the accuracy of the OTU generation step (test_otu.py)
2. For the resolution and accuracy of the classification of each source organism (test_classification.py) 

* Version
0.1

### How do I get set up? ###

One can download a working and fully installed version of decard in a docker container at:
https://hub.docker.com/r/golob/decard/
At that site is a tutorial to guide you through using decard. 

### Who do I talk to? ###

* jgolob@fhcrc.org