#!/bin/bash

#SBATCH -J create_genome_index
#SBATCH --partition comp06
#SBATCH -o Trimmomatic_%j.txt
#SBATCH -e Trimmomatic_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sdpotdar@uark.edu  
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=12
#SBATCH --time=01:00:00



### Convert genome gff annotation file to gtf using Cufflinks gffread ###
module load cufflinks

/share/apps/bioinformatics/cufflinks/cufflinks-2.2.1.Linux_x86_64/gffread \
/home/sdpotdar/H_melpomene_genome/Heliconius_melpomene_melpomene_Hmel2.5.gff3 \
-T \
-o /home/sdpotdar/H_melpomene_genome/Heliconius_melpomene_melpomene_Hmel2.5.gtf


### Create genome index ###
/home/sdpotdar/Bioinformatics_programs/STAR-2.7.1a/bin/Linux_x86_64_static/STAR \
--runThreadN 12 \
--runMode genomeGenerate \
--genomeDir /home/sdpotdar/H_melpomene_genome \
--genomeFastaFiles /home/sdpotdar/H_melpomene_genome/Heliconius_melpomene_melpomene_Hmel2.5.scaffolds.fa \
--sjdbGTFfile /home/sdpotdar/H_melpomene_genome/Heliconius_melpomene_melpomene_Hmel2.5.gtf \
--sjdbOverhang 49
