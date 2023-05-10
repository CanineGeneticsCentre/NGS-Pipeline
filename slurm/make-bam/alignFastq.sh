#!/usr/bin/env bash

#! sbatch directives begin here ###############################
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --time 04:00:00
#SBATCH --mail-type=BEGIN,END,FAIL,INVALID_DEPEND
#SBATCH -p skylake

#SBATCH -o logs/alignFastq_%A-%a.out

. /etc/profile.d/modules.sh                 # Leave this line (enables the module command)
module purge                                # Removes all modules still loaded
module load rhel7/default-peta4             # REQUIRED - loads the basic environment


SAMPLE=$1
LANE=$SLURM_ARRAY_TASK_ID
DIR=lane${LANE}

source ${SAMPLE}.config
mkdir $DIR

module load ${BWA}
module load ${SAMTOOLS}

SAMPLE=$1
LANE=$SLURM_ARRAY_TASK_ID
DIR=lane${LANE}

bwa mem -K 100000000 -p -v 3 -t 16 -Y ${FASTA}/${GENOME}.fasta ${DIR}/${SAMPLE}.L${LANE}.fastq.gz | samtools view -h -b -o ${DIR}/${SAMPLE}.L${LANE}.aligned.bam
