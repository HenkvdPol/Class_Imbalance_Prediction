#!/bin/bash
#SBATCH --job-name=sim
#SBATCH --output=logs/sim_%A_%a.out
#SBATCH --mail-user=h.van_der_pol@lumc.nl
#SBATCH --mail-type=ALL
#SBATCH --partition=gpu-l4-24g,gpu-2080ti-11g
#SBATCH --gres=gpu:1
#SBATCH --time=6-23:30:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --array=1-40

module purge
module load ALICE/default
module load R/4.4.2-gfbf-2024a
module load CMake/3.29.3-GCCcore-13.3.0

echo "========================================="
echo "Job array ID:  $SLURM_ARRAY_JOB_ID"
echo "Task index:    $SLURM_ARRAY_TASK_ID"
echo "Parameters:    $PARAMS"
echo "========================================="

PARAMS=$(sed -n "${SLURM_ARRAY_TASK_ID}p" parameters1.txt)
Rscript simulate.R $PARAMS
