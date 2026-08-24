#!/bin/bash
#SBATCH --job-name=boot_nomrdag
#SBATCH --time=08:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --array=1-4
#SBATCH --output=logs/nomrdag_idx%a_%A.out
#SBATCH --error=logs/nomrdag_idx%a_%A.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=peng.1276@osu.edu

module load r/4.5.1-5zezbzn

cd /home/peng.1276/matthew/MR_rr/sim_ArBr_bias_cluster
mkdir -p logs results

Rscript scripts/run_bootstrap_one.R ${SLURM_ARRAY_TASK_ID} no_mrdag