#!/usr/bin/env bash

#! sbatch directives begin here ###############################
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time 00:20:00
#SBATCH --mail-type=END,FAIL,INVALID_DEPEND
#SBATCH -p skylake

#SBATCH -o logs/qc-Yield_%A-%a.out

. /etc/profile.d/modules.sh                 # Leave this line (enables the module command)
module purge                                # Removes all modules still loaded
module load rhel7/default-peta4             # REQUIRED - loads the basic environment


SAMPLE=$1
LANE=$SLURM_ARRAY_TASK_ID
DIR=lane${LANE}
source ${SAMPLE}.config

module load ${GATK}

gatk CollectQualityYieldMetrics \
	--INPUT ${DIR}/${SAMPLE}.L${LANE}.unaligned.bam \
	--OUTPUT metrics/${SAMPLE}.L${LANE}.unmapped.quality_yield_metrics
