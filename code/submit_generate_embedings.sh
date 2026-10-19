#!/bin/bash
# Job Name
#PBS -N TOX_DMSO_UMAPs
#PBS -P je52
#PBS -q hugemem
#
# Resources
#PBS -l wd
#PBS -l ncpus=16
#PBS -l mem=1400GB
#PBS -l jobfs=500MB
#PBS -l walltime=5:00:00
#PBS -l storage=gdata/ei56+scratch/ei56+gdata/je52+scratch/je52+gdata/io72+gdata/wy83+scratch/io72+scratch/wy83
#PBS -j oe

#PBS -m ae
#PBS -M xochitl.happy@gmail.com

# For 2 files request 160 GB
# For 16 files request 1.3 TB
source ~/.bashrc

conda activate /g/data/je52/xd5809/envs/scAnalysis

python3 /g/data/wy83/cardiotoxicity/analysis/anndata/full_data_analysis/dim_reduction/code/generate_embedings_v2.py
