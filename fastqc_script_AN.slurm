#!/bin/bash

#SBATCH -J FastQC
#SBATCH --partition comp06
#SBATCH -o Fastqc_%j.txt
#SBATCH -e Fastqc_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sdpotdar@uark.edu  
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --time=6:00:00

export OMP_NUM_THREADS=32

### load fastqc module
module load java/sunjdk_1.8.0 fastqc

### Change directories to scratch space
cd /scratch/${SLURM_JOB_ID}

fastqc /scrfs/storage/sdpotdar/backup_data/AN_RAW_sequencing/EW-SP-41S-SP-01_S1_R1_001.fastq.gz \
fastqc /scrfs/storage/sdpotdar/backup_data/AN_RAW_sequencing/EW-SP-41S-SP-02_S2_R1_001.fastq.gz \
fastqc /scrfs/storage/sdpotdar/backup_data/AN_RAW_sequencing/EW-SP-41S-SP-03_S3_R1_001.fastq.gz \
fastqc /scrfs/storage/sdpotdar/backup_data/AN_RAW_sequencing/EW-SP-41S-SP-04_S4_R1_001.fastq.gz \
fastqc /scrfs/storage/sdpotdar/backup_data/AN_RAW_sequencing/EW-SP-41S-SP-05_S5_R1_001.fastq.gz \
fastqc /scrfs/storage/sdpotdar/backup_data/AN_RAW_sequencing/EW-SP-41S-SP-06_S6_R1_001.fastq.gz \
fastqc /scrfs/storage/sdpotdar/backup_data/AN_RAW_sequencing/EW-SP-41S-SP-07_S7_R1_001.fastq.gz \
fastqc /scrfs/storage/sdpotdar/backup_data/AN_RAW_sequencing/EW-SP-41S-SP-08_S8_R1_001.fastq.gz \
fastqc /scrfs/storage/sdpotdar/backup_data/AN_RAW_sequencing/EW-SP-41S-SP-09_S9_R1_001.fastq.gz \
fastqc /scrfs/storage/sdpotdar/backup_data/AN_RAW_sequencing/EW-SP-41S-SP-10_S10_R1_001.fastq.gz \
fastqc /scrfs/storage/sdpotdar/backup_data/AN_RAW_sequencing/EW-SP-41S-SP-11_S11_R1_001.fastq.gz \
fastqc /scrfs/storage/sdpotdar/backup_data/AN_RAW_sequencing/EW-SP-41S-SP-12_S12_R1_001.fastq.gz \
fastqc /scrfs/storage/sdpotdar/backup_data/AN_RAW_sequencing/EW-SP-41S-SP-13_S13_R1_001.fastq.gz \
fastqc /scrfs/storage/sdpotdar/backup_data/AN_RAW_sequencing/EW-SP-41S-SP-14_S14_R1_001.fastq.gz \
fastqc /scrfs/storage/sdpotdar/backup_data/AN_RAW_sequencing/EW-SP-41S-SP-15_S15_R1_001.fastq.gz \
fastqc /scrfs/storage/sdpotdar/backup_data/AN_RAW_sequencing/EW-SP-41S-SP-16_S16_R1_001.fastq.gz \
fastqc /scrfs/storage/sdpotdar/backup_data/AN_RAW_sequencing/EW-SP-41S-SP-17_S17_R1_001.fastq.gz \
fastqc /scrfs/storage/sdpotdar/backup_data/AN_RAW_sequencing/EW-SP-41S-SP-18_S18_R1_001.fastq.gz \
fastqc /scrfs/storage/sdpotdar/backup_data/AN_RAW_sequencing/EW-SP-41S-SP-19_S19_R1_001.fastq.gz \
fastqc /scrfs/storage/sdpotdar/backup_data/AN_RAW_sequencing/EW-SP-41S-SP-20_S20_R1_001.fastq.gz \
fastqc /scrfs/storage/sdpotdar/backup_data/AN_RAW_sequencing/EW-SP-41S-SP-21_S21_R1_001.fastq.gz \
fastqc /scrfs/storage/sdpotdar/backup_data/AN_RAW_sequencing/EW-SP-41S-SP-22_S22_R1_001.fastq.gz \
fastqc /scrfs/storage/sdpotdar/backup_data/AN_RAW_sequencing/EW-SP-41S-SP-23_S23_R1_001.fastq.gz \
fastqc /scrfs/storage/sdpotdar/backup_data/AN_RAW_sequencing/EW-SP-41S-SP-24_S24_R1_001.fastq.gz \
fastqc /scrfs/storage/sdpotdar/backup_data/AN_RAW_sequencing/EW-SP-41S-SP-25_S25_R1_001.fastq.gz \
fastqc /scrfs/storage/sdpotdar/backup_data/AN_RAW_sequencing/EW-SP-41S-SP-26_S26_R1_001.fastq.gz \
fastqc /scrfs/storage/sdpotdar/backup_data/AN_RAW_sequencing/EW-SP-41S-SP-27_S27_R1_001.fastq.gz \
fastqc /scrfs/storage/sdpotdar/backup_data/AN_RAW_sequencing/EW-SP-41S-SP-28_S28_R1_001.fastq.gz \
fastqc /scrfs/storage/sdpotdar/backup_data/AN_RAW_sequencing/EW-SP-41S-SP-29_S29_R1_001.fastq.gz \
fastqc /scrfs/storage/sdpotdar/backup_data/AN_RAW_sequencing/EW-SP-41S-SP-30_S30_R1_001.fastq.gz \
fastqc /scrfs/storage/sdpotdar/backup_data/AN_RAW_sequencing/EW-SP-41S-SP-31_S31_R1_001.fastq.gz \
fastqc /scrfs/storage/sdpotdar/backup_data/AN_RAW_sequencing/EW-SP-41S-SP-32_S32_R1_001.fastq.gz \
fastqc /scrfs/storage/sdpotdar/backup_data/AN_RAW_sequencing/EW-SP-41S-SP-33_S33_R1_001.fastq.gz \
fastqc /scrfs/storage/sdpotdar/backup_data/AN_RAW_sequencing/EW-SP-41S-SP-34_S34_R1_001.fastq.gz \
fastqc /scrfs/storage/sdpotdar/backup_data/AN_RAW_sequencing/EW-SP-41S-SP-35_S35_R1_001.fastq.gz \
fastqc /scrfs/storage/sdpotdar/backup_data/AN_RAW_sequencing/EW-SP-41S-SP-36_S36_R1_001.fastq.gz \
fastqc /scrfs/storage/sdpotdar/backup_data/AN_RAW_sequencing/EW-SP-41S-SP-37_S37_R1_001.fastq.gz \
fastqc /scrfs/storage/sdpotdar/backup_data/AN_RAW_sequencing/EW-SP-41S-SP-38_S38_R1_001.fastq.gz \
fastqc /scrfs/storage/sdpotdar/backup_data/AN_RAW_sequencing/EW-SP-41S-SP-39_S39_R1_001.fastq.gz \
fastqc /scrfs/storage/sdpotdar/backup_data/AN_RAW_sequencing/EW-SP-41S-SP-40_S40_R1_001.fastq.gz \
fastqc /scrfs/storage/sdpotdar/backup_data/AN_RAW_sequencing/EW-SP-41S-SP-41_S41_R1_001.fastq.gz \

mv *_fastqc.html /scrfs/storage/sdpotdar/backup_data/AN_raw_fastqc
mv *_fastqc.zip /scrfs/storage/sdpotdar/backup_data/AN_raw_fastqc
