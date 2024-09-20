#!/bin/bash
#SBATCH -J 02_1.munge
#SBATCH -p cpu
#SBATCH -o ./%x_%j.out
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=14-00:00:00
#SBATCH --mail-user=sh.austin.park@gmail.com
#SBATCH --mail-type=END,FAIL

module purge

CONDA_PATH=/data1/software/anaconda3
ENV_NAME=ldsc
ENV_PATH=$CONDA_PATH/envs/$ENV_NAME
source $CONDA_PATH/bin/activate $ENV_PATH

conda activate ldsc
### ### ### ### ### ### ### ### ### ### ### ### ###

metadata_prefix="metadata_for_external_trait_rg"

### ### ### ### ### ### ### ### ### ### ### ### ###
src_munge=/data1/sanghyeon/wonlab_contribute/combined/software/ldsc/munge_sumstats.py
SNPlist=/data1/sanghyeon/wonlab_contribute/combined/software/ldsc/data/w_hm3.snplist

f_metadata="${metadata_prefix}.for_munging.csv"
dir_out="./munged"; mkdir -p ${dir_out} 

tail -n +2 "${f_metadata}" | while IFS=',' read -r trait f_gwas; do
    if [ ! -f "${dir_out}/${trait}.sumstats.gz" ]; then
        /data1/software/anaconda3/envs/ldsc/bin/python2 ${src_munge} \
            --sumstats ${f_gwas} \
            --p P \
            --out "${dir_out}/${trait}" \
            --merge-alleles ${SNPlist} \
            --chunksize 500000
    fi
done

