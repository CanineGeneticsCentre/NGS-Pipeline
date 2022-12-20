#!/usr/bin/env bash

#SBATCH -A MELLERSH-SL3-CPU
#SBATCH -o /home/%u/hpc-work/logs/job-%j.out
#SBATCH -p icelake
#SBATCH -t 01:00:00
#SBATCH --mail-type=FAIL,INVALID_DEPEND,END


FROM=$1
TO=$2

echo rsync --progress -auvh ${FROM} ${TO}