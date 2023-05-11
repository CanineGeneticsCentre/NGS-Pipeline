#!/usr/bin/env bash

#! sbatch directives begin here ###############################
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time 01:00:00
#SBATCH --mail-type=END,FAIL,INVALID_DEPEND
#SBATCH -p cclake

#SBATCH -o logs/qc-Lane_%A-%a.out

. /etc/profile.d/modules.sh                 # Leave this line (enables the module command)
module purge
module load rhel7/default-ccl


SAMPLE=$1
LANE=$SLURM_ARRAY_TASK_ID
DIR=lane${LANE}
source ${SAMPLE}.config

module load $GATK

gatk CollectMultipleMetrics \
	--INPUT ${DIR}/${SAMPLE}.L${LANE}.merged.bam \
	--OUTPUT metrics/${SAMPLE}.L${LANE}.readgroup \
	--ASSUME_SORTED true \
	--PROGRAM null \
	--PROGRAM CollectBaseDistributionByCycle \
	--PROGRAM CollectInsertSizeMetrics \
	--PROGRAM MeanQualityByCycle \
	--PROGRAM QualityScoreDistribution \
	--METRIC_ACCUMULATION_LEVEL null \
	--METRIC_ACCUMULATION_LEVEL ALL_READS
