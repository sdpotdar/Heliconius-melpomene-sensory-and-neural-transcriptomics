#!/bin/bash

### Extract gene sequences only from gff3 file
awk '$3 == "gene" {print $0}'  ~/Documents/uark/phd/projects/melpomene_rnaseq/bioinformatics/H_melpomene_genome/Hmel2.5.gff3 > ~/Documents/uark/phd/projects/melpomene_rnaseq/bioinformatics/H_melpomene_genome/Hmel2.5_genes_only.gff3

### Switch column 9 (containing gene ID) with column 3 for use with BedTools getfasta
awk 'OFS="\t" {print $1,$2,$9,$4,$5,$6,$7,$8,$3}'  ~/Documents/uark/phd/projects/melpomene_rnaseq/bioinformatics/H_melpomene_genome/Hmel2.5_genes_only.gff3 > ~/Documents/uark/phd/projects/melpomene_rnaseq/bioinformatics/H_melpomene_genome/Hmel2.5_genes_only_column_reorder.gff3

### Remove ID= from column 3, leaving just gene ID
awk '{gsub("ID=","");print}' ~/Documents/uark/phd/projects/melpomene_rnaseq/bioinformatics/H_melpomene_genome/Hmel2.5_genes_only_column_reorder.gff3 > ~/Documents/uark/phd/projects/melpomene_rnaseq/bioinformatics/H_melpomene_genome/Hmel2.5_genes_only_column_reorder_IDs.gff3

### Make fasta of Heliconius genes with Bedtools getfasta
~/Downloads/bedtools2/bin/bedtools \
getfasta \
-fi ~/Documents/uark/phd/projects/melpomene_rnaseq/bioinformatics/H_melpomene_genome/Heliconius_melpomene_melpomene_Hmel2.5.cds.fa \
-bed ~/Documents/uark/phd/projects/melpomene_rnaseq/bioinformatics/H_melpomene_genome/Hmel2.5_genes_only_column_reorder_IDs.gff3 \
-fo ~/Documents/uark/phd/projects/melpomene_rnaseq/bioinformatics/H_melpomene_genome/Hmel_2.5_genes_only.fa \
-name