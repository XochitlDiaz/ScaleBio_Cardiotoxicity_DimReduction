#!/bin/bash
# Job Name
#PBS -N full_data_prep
#PBS -P wy83
#PBS -q normal
# ### rffe PBS -q hugemem
#
# Resources
#PBS -l wd
#PBS -l ncpus=16
#PBS -l mem=160GB
#PBS -l jobfs=500MB
#PBS -l walltime=48:00:00
#PBS -l storage=gdata/ei56+scratch/ei56+gdata/je52+scratch/je52+gdata/io72+gdata/wy83+scratch/io72+scratch/wy83
#PBS -j oe

#PBS -m ae
#PBS -M xochitl.happy@gmail.com

# For 2 files request 160 GB
# For 16 files request 1.3 TB
source ~/.bashrc

conda activate /g/data/je52/xd5809/envs/scAnalysis
cd /g/data/wy83/cardiotoxicity/analysis/anndata/full_data_analysis/dim_reduction/code
python3 -m pipeline.run_pipeline
