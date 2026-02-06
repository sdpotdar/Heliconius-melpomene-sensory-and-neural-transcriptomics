#!/bin/bash
#SBATCH -J HE_CONCAT_DIAMOND_BLASTX
#SBATCH --partition comp06
#SBATCH -o header_edited_concat_part-1_diamond_blastx_%j.txt
#SBATCH -e header_edited_concat_part-1_diamond_blastx_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --time=06:00:00

#Load Diamond module
module load diamond/2.1.7
cd diamond_blast

echo "Blasting sequence $f"
diamond blastx \
        --tmpdir /scratch/$SLURM_JOB_ID \
        --threads 32 \
        --query /scrfs/storage/sdpotdar/backup_data/Hmel_2.5_genes_only_concatenated_header_edited/Hmel_2.5_genes_only_header_edited_part-1.fa \
        --out /scrfs/storage/sdpotdar/backup_data/Concat_blastx_header_edited/Hmel_2.5_genes_header_edited_part-1.xml \
        --db /scrfs/storage/sdpotdar/fast/nr.dmnd \
        --evalue 1e-3 \
        --outfmt 5 \
        --max-target-seqs 10

echo "Finished sequence $f"

rsync -a diamond_blast /scrfs/storage/sdpotdar/backup_data/
### Each part has a thousand seqs. Part-1 done ###
### 20097 total sequences ###
