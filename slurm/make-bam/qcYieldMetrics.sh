#!/usr/bin/env bash

#! sbatch directives begin here ###############################
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time 00:20:00
#SBATCH --mail-type=FAIL
#SBATCH -p cclake

#SBATCH -o logs/qc-Yield_%A-%a.out

. /etc/profile.d/modules.sh                 # Leave this line (enables the module command)
module purge
module load rhel7/default-ccl


SAMPLE=$1
LANE=$SLURM_ARRAY_TASK_ID
DIR=lane${LANE}
source ${SAMPLE}.config

module load ${GATK}

gatk CollectQualityYieldMetrics \
	--INPUT ${DIR}/${SAMPLE}.L${LANE}.unaligned.bam \
	--OUTPUT metrics/${SAMPLE}.L${LANE}.unmapped.quality_yield_metrics
