#!/bin/sh
#SBATCH --ntasks=4
#SBATCH --nodes=1
#SBATCH --mem=30gb
#SBATCH --time=1:00:00
#SBATCH --partition=batch
#SBATCH --licenses=common
#SBATCH --job-name=mediation_demo
#SBATCH --mail-user=zyang35@huskers.unl.edu
#SBATCH --mail-type=ALL
#SBATCH --error=/common/jyanglab/zhikaiyang/projects/high_dim_mediation/slurm-log/mediation_demo.err
#SBATCH --output=/common/jyanglab/zhikaiyang/projects/high_dim_mediation/slurm-log/mediation_demo.out
module load R/3.5
Rscript /common/jyanglab/zhikaiyang/projects/high_dim_mediation/workflow/1_mediation_demo.R
