#!/usr/bin/env bash

#SBATCH -A MELLERSH-SL3-CPU
#SBATCH -o /home/%u/hpc-work/logs/job-%j.out
#SBATCH -p icelake
#SBATCH -t 04:00:00
#SBATCH --mail-type=FAIL,INVALID_DEPEND,END


source $CONDA_PREFIX/etc/profile.d/conda.sh  # Always add this command to your scripts
conda activate ENA

ascp -QT -l300M -L- *.fq.gz Webin-47111@webin.ebi.ac.uk:.