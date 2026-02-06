#!/bin/bash

### Split fasta into single seqs ###
perl ~/Downloads/fasta-splitter-0.2.6/fasta-splitter.pl \
--part-size 1000 \
--measure count \
--line-length 0 \
--nopad \
--out-dir ~/Documents/uark/phd/projects/melpomene_rnaseq/bioinformatics/H_melpomene_genome/Hmel_2.5_genes_only.fa_split/thousand_seqs \
~/Documents/uark/phd/projects/melpomene_rnaseq/bioinformatics/H_melpomene_genome/Hmel_2.5_genes_only.fa