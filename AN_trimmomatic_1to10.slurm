#!/bin/bash
#SBATCH -J trimmomatic_hmel
#SBATCH --partition comp06
#SBATCH -o Trimmomatic_Hmel_%j.txt
#SBATCH -e Trimmomatic_Hmel_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sdpotdar@uark.edu  
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --time=06:00:00
export OMP_NUM_THREADS=32

### Load required modules
module load java/sunjdk_1.8.0

### Change directories to scratch space
cd /scratch/${SLURM_JOB_ID}

cp /scrfs/storage/sdpotdar/backup_data/AN_RAW_sequencing/EW-SP-41S-SP-01_S1_R1_001.fastq.gz .

### Remove adapters using Trimmomatic
echo Removing Illumina adapters with Trimmomatic...
java -jar /share/apps/bioinformatics/trimmomatic/Trimmomatic-0.36/trimmomatic-0.36.jar \
SE -phred33 \
EW-SP-41S-SP-01_S1_R1_001.fastq.gz \
EW-SP-41S-SP-01_S1_R1_001.fastq.gz_adapt-trimmed.fastq.gz \
ILLUMINACLIP:/share/apps/bioinformatics/trimmomatic/Trimmomatic-0.36/adapters/TruSeq3-SE.fa:2:30:10 MINLEN:100

cp /scrfs/storage/sdpotdar/backup_data/AN_RAW_sequencing/EW-SP-41S-SP-02_S2_R1_001.fastq.gz .

### Remove adapters using Trimmomatic
echo Removing Illumina adapters with Trimmomatic...
java -jar /share/apps/bioinformatics/trimmomatic/Trimmomatic-0.36/trimmomatic-0.36.jar \
SE -phred33 \
EW-SP-41S-SP-02_S2_R1_001.fastq.gz \
EW-SP-41S-SP-02_S2_R1_001.fastq.gz_adapt-trimmed.fastq.gz \
ILLUMINACLIP:/share/apps/bioinformatics/trimmomatic/Trimmomatic-0.36/adapters/TruSeq3-SE.fa:2:30:10 MINLEN:100

cp /scrfs/storage/sdpotdar/backup_data/AN_RAW_sequencing/EW-SP-41S-SP-03_S3_R1_001.fastq.gz .

### Remove adapters using Trimmomatic
echo Removing Illumina adapters with Trimmomatic...
java -jar /share/apps/bioinformatics/trimmomatic/Trimmomatic-0.36/trimmomatic-0.36.jar \
SE -phred33 \
EW-SP-41S-SP-03_S3_R1_001.fastq.gz \
EW-SP-41S-SP-03_S3_R1_001.fastq.gz_adapt-trimmed.fastq.gz \
ILLUMINACLIP:/share/apps/bioinformatics/trimmomatic/Trimmomatic-0.36/adapters/TruSeq3-SE.fa:2:30:10 MINLEN:100

cp /scrfs/storage/sdpotdar/backup_data/AN_RAW_sequencing/EW-SP-41S-SP-04_S4_R1_001.fastq.gz .

### Remove adapters using Trimmomatic
echo Removing Illumina adapters with Trimmomatic...
java -jar /share/apps/bioinformatics/trimmomatic/Trimmomatic-0.36/trimmomatic-0.36.jar \
SE -phred33 \
EW-SP-41S-SP-04_S4_R1_001.fastq.gz \
EW-SP-41S-SP-04_S4_R1_001.fastq.gz_adapt-trimmed.fastq.gz \
ILLUMINACLIP:/share/apps/bioinformatics/trimmomatic/Trimmomatic-0.36/adapters/TruSeq3-SE.fa:2:30:10 MINLEN:100

cp /scrfs/storage/sdpotdar/backup_data/AN_RAW_sequencing/EW-SP-41S-SP-05_S5_R1_001.fastq.gz .

### Remove adapters using Trimmomatic
echo Removing Illumina adapters with Trimmomatic...
java -jar /share/apps/bioinformatics/trimmomatic/Trimmomatic-0.36/trimmomatic-0.36.jar \
SE -phred33 \
EW-SP-41S-SP-05_S5_R1_001.fastq.gz \
EW-SP-41S-SP-05_S5_R1_001.fastq.gz_adapt-trimmed.fastq.gz \
ILLUMINACLIP:/share/apps/bioinformatics/trimmomatic/Trimmomatic-0.36/adapters/TruSeq3-SE.fa:2:30:10 MINLEN:100

cp /scrfs/storage/sdpotdar/backup_data/AN_RAW_sequencing/EW-SP-41S-SP-06_S6_R1_001.fastq.gz .

### Remove adapters using Trimmomatic
echo Removing Illumina adapters with Trimmomatic...
java -jar /share/apps/bioinformatics/trimmomatic/Trimmomatic-0.36/trimmomatic-0.36.jar \
SE -phred33 \
EW-SP-41S-SP-06_S6_R1_001.fastq.gz \
EW-SP-41S-SP-06_S6_R1_001.fastq.gz_adapt-trimmed.fastq.gz \
ILLUMINACLIP:/share/apps/bioinformatics/trimmomatic/Trimmomatic-0.36/adapters/TruSeq3-SE.fa:2:30:10 MINLEN:100

cp /scrfs/storage/sdpotdar/backup_data/AN_RAW_sequencing/EW-SP-41S-SP-07_S7_R1_001.fastq.gz .

### Remove adapters using Trimmomatic
echo Removing Illumina adapters with Trimmomatic...
java -jar /share/apps/bioinformatics/trimmomatic/Trimmomatic-0.36/trimmomatic-0.36.jar \
SE -phred33 \
EW-SP-41S-SP-07_S7_R1_001.fastq.gz \
EW-SP-41S-SP-07_S7_R1_001.fastq.gz_adapt-trimmed.fastq.gz \
ILLUMINACLIP:/share/apps/bioinformatics/trimmomatic/Trimmomatic-0.36/adapters/TruSeq3-SE.fa:2:30:10 MINLEN:100

cp /scrfs/storage/sdpotdar/backup_data/AN_RAW_sequencing/EW-SP-41S-SP-08_S8_R1_001.fastq.gz .

### Remove adapters using Trimmomatic
echo Removing Illumina adapters with Trimmomatic...
java -jar /share/apps/bioinformatics/trimmomatic/Trimmomatic-0.36/trimmomatic-0.36.jar \
SE -phred33 \
EW-SP-41S-SP-08_S8_R1_001.fastq.gz \
EW-SP-41S-SP-08_S8_R1_001.fastq.gz_adapt-trimmed.fastq.gz \
ILLUMINACLIP:/share/apps/bioinformatics/trimmomatic/Trimmomatic-0.36/adapters/TruSeq3-SE.fa:2:30:10 MINLEN:100

cp /scrfs/storage/sdpotdar/backup_data/AN_RAW_sequencing/EW-SP-41S-SP-09_S9_R1_001.fastq.gz .

### Remove adapters using Trimmomatic
echo Removing Illumina adapters with Trimmomatic...
java -jar /share/apps/bioinformatics/trimmomatic/Trimmomatic-0.36/trimmomatic-0.36.jar \
SE -phred33 \
EW-SP-41S-SP-09_S9_R1_001.fastq.gz \
EW-SP-41S-SP-09_S9_R1_001.fastq.gz_adapt-trimmed.fastq.gz \
ILLUMINACLIP:/share/apps/bioinformatics/trimmomatic/Trimmomatic-0.36/adapters/TruSeq3-SE.fa:2:30:10 MINLEN:100

cp /scrfs/storage/sdpotdar/backup_data/AN_RAW_sequencing/EW-SP-41S-SP-10_S10_R1_001.fastq.gz .

### Remove adapters using Trimmomatic
echo Removing Illumina adapters with Trimmomatic...
java -jar /share/apps/bioinformatics/trimmomatic/Trimmomatic-0.36/trimmomatic-0.36.jar \
SE -phred33 \
EW-SP-41S-SP-10_S10_R1_001.fastq.gz \
EW-SP-41S-SP-10_S10_R1_001.fastq.gz_adapt-trimmed.fastq.gz \
ILLUMINACLIP:/share/apps/bioinformatics/trimmomatic/Trimmomatic-0.36/adapters/TruSeq3-SE.fa:2:30:10 MINLEN:100
