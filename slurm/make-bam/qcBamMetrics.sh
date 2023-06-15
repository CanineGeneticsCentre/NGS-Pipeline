#!/usr/bin/env bash

#! RUN : sbatch bamMetrics.sh <SAMPLE> <REF>

#! sbatch directives begin here ###############################
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time 12:00:00
#SBATCH --mail-type=FAIL,END
#SBATCH -p cclake

#SBATCH -o logs/bamMetrics-%j.out

. /etc/profile.d/modules.sh                 # Leave this line (enables the module command)
module purge
module load rhel7/default-ccl


SAMPLE=$1
REF=$2
source ${SAMPLE}.config

module load ${SAMTOOLS}
module load $PICARD

BAM_FILE=${SAMPLE}-${REF}.bam


samtools flagstat ${BAM_FILE} > metrics/${SAMPLE}.flagstat
picard_latest CollectWgsMetrics INPUT=${BAM_FILE} O=metrics/${SAMPLE}.wgs.txt REFERENCE_SEQUENCE=${FASTA}/${GENOME}.fasta STOP_AFTER=100000000 VALIDATION_STRINGENCY=LENIENT
picard_latest CollectInsertSizeMetrics INPUT=${BAM_FILE} O=metrics/${SAMPLE}.insert_size.txt H=metrics/${SAMPLE}.insert_size_histogram.pdf VALIDATION_STRINGENCY=LENIENT

echo ${SAMPLE}
echo "..."

grep 'in total' metrics/${SAMPLE}.flagstat
grep '%' metrics/${SAMPLE}.flagstat | grep -v 'singleton'
echo "..."

head -8 metrics/${SAMPLE}.insert_size.txt | tail -2 | cut -f 6,8
echo "..."

head -8 metrics/${SAMPLE}.wgs.txt | tail -2 | cut -f 2,3,15,16,18,20-22
echo "..."
