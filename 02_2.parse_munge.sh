#!/bin/bash
#SBATCH -J 02_2.parse_munge
#SBATCH -p cpu
#SBATCH -o ./%x_%j.out
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --time=14-00:00:00
#SBATCH --mail-user=sh.austin.park@gmail.com 
#SBATCH --mail-type=END,FAIL 

module purge

CONDA_PATH=/data1/software/anaconda3 
source $CONDA_PATH/bin/activate /data1/software/anaconda3/envs/R_SP_4.4 

/data1/software/anaconda3/envs/R_SP_4.4/bin/Rscript 02_2.parse_munge.R