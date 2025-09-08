#!/bin/bash
#SBATCH -n 32
#SBATCH -N 1
#SBATCH --mem 15000
#SBATCH -o kallisto_%A_%a.out
#SBATCH -e kallisto_%A_%a.err
#SBATCH -J kallisto_arr
#SBATCH -t 120
#SBATCH --array=1-3

module load gcc openmpi kallisto

TRANS=/gpfs/scratch/finer03
KALLISTO_DIR=/gpfs/scratch/finer03/kallisto

dataset1=HAP1_${SLURM_ARRAY_TASK_ID}
dataset2=eHAP1_${SLURM_ARRAY_TASK_ID}

kallisto quant -i /gpfs/home/finer03/homo_sapiens/transcriptome.idx -o $KALLISTO_DIR/${dataset1}.kallisto_out -b 100 -t 32 $TRANS/${dataset1}_1.fastq.gz $TRANS/${dataset1}_2.fastq.gz

kallisto quant -i /gpfs/home/finer03/homo_sapiens/transcriptome.idx -o $KALLISTO_DIR/${dataset2}.kallisto_out -b 100 -t 32 $TRANS/${dataset2}_1.fastq.gz $TRANS/${dataset2}_2.fastq.gz

