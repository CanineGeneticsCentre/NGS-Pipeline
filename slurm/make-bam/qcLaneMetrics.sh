#!/usr/bin/env bash

#! sbatch directives begin here ###############################
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time 01:00:00
#SBATCH --mail-type=FAIL
#SBATCH -p cclake

#SBATCH -o logs/qc-Lane_%A-%a.out

. /etc/profile.d/modules.sh                # Leave this line (enables the module command)
module purge                               # Removes all modules still loaded
module load rhel8/default-ccl              # REQUIRED - loads the basic environment


SAMPLE=$1
DIR=rg${SLURM_ARRAY_TASK_ID}
source ${SAMPLE}.config

module load $GATK

gatk CollectMultipleMetrics \
	--INPUT ${DIR}/${SAMPLE}.RG${SLURM_ARRAY_TASK_ID}.merged.bam \
	--OUTPUT metrics/${SAMPLE}.RG${SLURM_ARRAY_TASK_ID}.readgroup \
	--ASSUME_SORTED true \
	--PROGRAM null \
	--PROGRAM CollectBaseDistributionByCycle \
	--PROGRAM CollectInsertSizeMetrics \
	--PROGRAM MeanQualityByCycle \
	--PROGRAM QualityScoreDistribution \
	--METRIC_ACCUMULATION_LEVEL null \
	--METRIC_ACCUMULATION_LEVEL ALL_READS
