#!/usr/bin/env bash

#SBATCH -A MELLERSH-SL3-CPU
#SBATCH -p cclake
#SBATCH -t 08:00:00
#SBATCH --mail-type=FAIL,INVALID_DEPEND,END
##SBATCH --mail-type=ALL
##SBATCH -o logs/ena-%j.out
#SBATCH -o /home/%u/hpc-work/logs/ena-%A_%a.out


source $CONDA_PREFIX/etc/profile.d/conda.sh  # Always add this command to your scripts
conda activate ENA

SAMPLE=`head -${SLURM_ARRAY_TASK_ID} ${SAMPLE_LIST} | tail -1`

DIR=$1
cd ${SAMPLE}/${DIR}

ascp -QT -l300M -k2 -L- * Webin-47111@webin.ebi.ac.uk:.
#ascp -QT -l300M -L- * Webin-47111@webin.ebi.ac.uk:.

#cd ../; rm -rf $DIR; cd ../
