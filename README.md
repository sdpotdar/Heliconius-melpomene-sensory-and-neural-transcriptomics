# Heliconius-melpomene-sensory-and-neural-transcriptomics
This repository contains the codes (and the parameters within them) to identify the neurotranscriptomic signatures associated with natural variation in mate preference learning in two subspecies of Heliconius melpomene. The example codes are for the antennae (AN) tissues, but the same codes and parameters were used for the eyes (EY) and the brain (BR) tissues.

The sequence of the analyses were as follows:

A. Bioinformatics to get read counts of genes from raw reads
1. Check quality of raw reads using fastqc: fastqc_script_AN.slurm
2. Remove adapters from raw reads using Trimmomatic v0.38: AN_trimmomatic_1to10.slurm
3. Create STAR genome index for H. melpomene v2.5 genome: create_STAR_genome_index_Hmel.slurm
4. Align trimmed reads to the genome using STAR: AN_SP01to10_align_reads_to_genome.slurm
 
B. Build DESeq datasets to analyse differential gene expression between groups and treatments
1. Use DESeq2 to identify all the expressed genes in each tissue: Hmel_Antenna_AEG.R
2. Use DESeq2 to analyse Differential Gene Expression between various groups and treatments: Hmel_DEseq_AN.R
3. Add annotations to DEG/AEG files: DEG_annotation.R
4. Add interproscan data and gene ontology IDs to DEG/AEG files: merge_interproscan_genelist.R
5. Identify chromosome specific DEGs: DEG_Chr_location.R
6. Plot heatmaps for all DEGs: Pub_AN_heatmaps.R
7. Plot specific gene DEG barplots and PCA plots: plot_counts_figures.R
8. Make venn diagram plots for unique and shared DEGs between treatments: venn_plots.R

C. Build functional annotation files for Blast2GO and perform Gene Ontology Enrichment Analyses
1. Create genome fasta files for genes only: make_genome_fasta_with_bedtools_hmelpomene.sh
2. Split Hmel genes into 21 fasta file, each containing 1000 sequences: split_Hmel_genes_fasta.sh
3. Make Diamond database: make_diamond_db.slurm
4. Blast sequences against diamond database: diamond_blast_with_cat_header_edited_part-1.slurm

The following .xml files can now be uploaded to Blast2GO for GUI based GO enrichment and functional analyses.

D. Weighted Gene Co-expression Network Analysis (WGCNA)
1. Perform WGCNA after running DESeq2: WGCNA_AN.R
2. Extract important modules after running WGCNA: WGCNA_extract_data_for_important_modules--ANTENNA-SIGNED.R
