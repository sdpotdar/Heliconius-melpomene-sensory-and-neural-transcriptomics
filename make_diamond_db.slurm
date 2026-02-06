#!/bin/bash
#SBATCH -J DIAMOND_DB
#SBATCH --partition comp06
#SBATCH -o make_diamond_db_%j.txt
#SBATCH -e make_diamond_db_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=06:00:00

module load gcc-11.2.1/SKYLAKEX/diamond/2.0.15

#Download blast fasta nr database for usage in diamond
mkdir /scrfs/storage/sdpotdar/backup_data/diamond_db
wget 'ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz'

#Run makedb for diamond using the downloaded nr.gz from ncbi
diamond makedb --in nr.gz -d nr
echo "finished making nr.gz database"
