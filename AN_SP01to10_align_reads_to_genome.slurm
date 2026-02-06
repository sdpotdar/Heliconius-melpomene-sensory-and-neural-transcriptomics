#!/bin/bash
#SBATCH -J align_reads_genome_STAR_Hmel
#SBATCH --partition comp06
#SBATCH -o alignment_%j.txt
#SBATCH -e alignment_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sdpotdar@uark.edu  
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --time=06:00:00

export OMP_NUM_THREADS=32
module purge
module load STAR

###Change directories
cd /scrfs/storage/sdpotdar/backup_data/AN_trimmed_fastq


### Align reads to genome ###
STAR --runThreadN 32 \
--genomeDir /home/sdpotdar/H_melpomene_genome \
--sjdbGTFfile /home/sdpotdar/H_melpomene_genome/Heliconius_melpomene_melpomene_Hmel2.5.gtf \
--sjdbOverhang 99 \
--quantMode GeneCounts \
--twopassMode Basic \
--outReadsUnmapped Fastx \
--readFilesIn /scrfs/storage/sdpotdar/backup_data/AN_trimmed_fastq/EW-SP-41S-SP-01_S1_R1_001.fastq.gz_adapt-trimmed.fastq.gz \
--outFileNamePrefix  SP-01 \
--readFilesCommand zcat

### Align reads to genome ###
STAR --runThreadN 32 \
--genomeDir /home/sdpotdar/H_melpomene_genome \
--sjdbGTFfile /home/sdpotdar/H_melpomene_genome/Heliconius_melpomene_melpomene_Hmel2.5.gtf \
--sjdbOverhang 99 \
--quantMode GeneCounts \
--twopassMode Basic \
--outReadsUnmapped Fastx \
--readFilesIn /scrfs/storage/sdpotdar/backup_data/AN_trimmed_fastq/EW-SP-41S-SP-02_S2_R1_001.fastq.gz_adapt-trimmed.fastq.gz \
--outFileNamePrefix  SP-02 \
--readFilesCommand zcat

### Align reads to genome ###
STAR --runThreadN 32 \
--genomeDir /home/sdpotdar/H_melpomene_genome \
--sjdbGTFfile /home/sdpotdar/H_melpomene_genome/Heliconius_melpomene_melpomene_Hmel2.5.gtf \
--sjdbOverhang 99 \
--quantMode GeneCounts \
--twopassMode Basic \
--outReadsUnmapped Fastx \
--readFilesIn /scrfs/storage/sdpotdar/backup_data/AN_trimmed_fastq/EW-SP-41S-SP-03_S3_R1_001.fastq.gz_adapt-trimmed.fastq.gz \
--outFileNamePrefix  SP-03 \
--readFilesCommand zcat

### Align reads to genome ###
STAR --runThreadN 32 \
--genomeDir /home/sdpotdar/H_melpomene_genome \
--sjdbGTFfile /home/sdpotdar/H_melpomene_genome/Heliconius_melpomene_melpomene_Hmel2.5.gtf \
--sjdbOverhang 99 \
--quantMode GeneCounts \
--twopassMode Basic \
--outReadsUnmapped Fastx \
--readFilesIn /scrfs/storage/sdpotdar/backup_data/AN_trimmed_fastq/EW-SP-41S-SP-04_S4_R1_001.fastq.gz_adapt-trimmed.fastq.gz \
--outFileNamePrefix  SP-04 \
--readFilesCommand zcat

### Align reads to genome ###
STAR --runThreadN 32 \
--genomeDir /home/sdpotdar/H_melpomene_genome \
--sjdbGTFfile /home/sdpotdar/H_melpomene_genome/Heliconius_melpomene_melpomene_Hmel2.5.gtf \
--sjdbOverhang 99 \
--quantMode GeneCounts \
--twopassMode Basic \
--outReadsUnmapped Fastx \
--readFilesIn /scrfs/storage/sdpotdar/backup_data/AN_trimmed_fastq/EW-SP-41S-SP-05_S5_R1_001.fastq.gz_adapt-trimmed.fastq.gz \
--outFileNamePrefix  SP-05 \
--readFilesCommand zcat

### Align reads to genome ###
STAR --runThreadN 32 \
--genomeDir /home/sdpotdar/H_melpomene_genome \
--sjdbGTFfile /home/sdpotdar/H_melpomene_genome/Heliconius_melpomene_melpomene_Hmel2.5.gtf \
--sjdbOverhang 99 \
--quantMode GeneCounts \
--twopassMode Basic \
--outReadsUnmapped Fastx \
--readFilesIn /scrfs/storage/sdpotdar/backup_data/AN_trimmed_fastq/EW-SP-41S-SP-06_S6_R1_001.fastq.gz_adapt-trimmed.fastq.gz \
--outFileNamePrefix  SP-06 \
--readFilesCommand zcat

### Align reads to genome ###
STAR --runThreadN 32 \
--genomeDir /home/sdpotdar/H_melpomene_genome \
--sjdbGTFfile /home/sdpotdar/H_melpomene_genome/Heliconius_melpomene_melpomene_Hmel2.5.gtf \
--sjdbOverhang 99 \
--quantMode GeneCounts \
--twopassMode Basic \
--outReadsUnmapped Fastx \
--readFilesIn /scrfs/storage/sdpotdar/backup_data/AN_trimmed_fastq/EW-SP-41S-SP-07_S7_R1_001.fastq.gz_adapt-trimmed.fastq.gz \
--outFileNamePrefix  SP-07 \
--readFilesCommand zcat

### Align reads to genome ###
STAR --runThreadN 32 \
--genomeDir /home/sdpotdar/H_melpomene_genome \
--sjdbGTFfile /home/sdpotdar/H_melpomene_genome/Heliconius_melpomene_melpomene_Hmel2.5.gtf \
--sjdbOverhang 99 \
--quantMode GeneCounts \
--twopassMode Basic \
--outReadsUnmapped Fastx \
--readFilesIn /scrfs/storage/sdpotdar/backup_data/AN_trimmed_fastq/EW-SP-41S-SP-08_S8_R1_001.fastq.gz_adapt-trimmed.fastq.gz \
--outFileNamePrefix  SP-08 \
--readFilesCommand zcat

### Align reads to genome ###
STAR --runThreadN 32 \
--genomeDir /home/sdpotdar/H_melpomene_genome \
--sjdbGTFfile /home/sdpotdar/H_melpomene_genome/Heliconius_melpomene_melpomene_Hmel2.5.gtf \
--sjdbOverhang 99 \
--quantMode GeneCounts \
--twopassMode Basic \
--outReadsUnmapped Fastx \
--readFilesIn /scrfs/storage/sdpotdar/backup_data/AN_trimmed_fastq/EW-SP-41S-SP-09_S9_R1_001.fastq.gz_adapt-trimmed.fastq.gz \
--outFileNamePrefix  SP-09 \
--readFilesCommand zcat

### Align reads to genome ###
STAR --runThreadN 32 \
--genomeDir /home/sdpotdar/H_melpomene_genome \
--sjdbGTFfile /home/sdpotdar/H_melpomene_genome/Heliconius_melpomene_melpomene_Hmel2.5.gtf \
--sjdbOverhang 99 \
--quantMode GeneCounts \
--twopassMode Basic \
--outReadsUnmapped Fastx \
--readFilesIn /scrfs/storage/sdpotdar/backup_data/AN_trimmed_fastq/EW-SP-41S-SP-10_S10_R1_001.fastq.gz_adapt-trimmed.fastq.gz \
--outFileNamePrefix  SP-10 \
--readFilesCommand zcat