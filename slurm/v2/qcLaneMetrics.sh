#!/usr/bin/env bash

#! sbatch directives begin here ###############################
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time 01:00:00
#SBATCH --mail-type=END,FAIL,INVALID_DEPEND
#SBATCH -p skylake

#SBATCH -o logs/qc-Lane_%j.out

. /etc/profile.d/modules.sh                 # Leave this line (enables the module command)
module purge                                # Removes all modules still loaded
module load rhel7/default-peta4             # REQUIRED - loads the basic environment

module load gatk-4.2.5.0-gcc-5.4.0-hzdcjga

LANE=$1
DIR=lane${LANE}

gatk CollectMultipleMetrics \
	--INPUT ${DIR}/ESD_37411.L${LANE}.merged.bam \
	--OUTPUT metrics/ESD_37411.L${LANE}.readgroup \
	--ASSUME_SORTED true \
	--PROGRAM null \
	--PROGRAM CollectBaseDistributionByCycle \
	--PROGRAM CollectInsertSizeMetrics \
	--PROGRAM MeanQualityByCycle \
	--PROGRAM QualityScoreDistribution \
	--METRIC_ACCUMULATION_LEVEL null \
	--METRIC_ACCUMULATION_LEVEL ALL_READS
