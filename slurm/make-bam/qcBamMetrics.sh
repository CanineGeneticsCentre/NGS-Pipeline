#!/usr/bin/env bash

#! RUN : sbatch bamMetrics.sh <SAMPLE> <REF>

#! sbatch directives begin here ###############################
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time 12:00:00
#SBATCH --mail-type=FAIL,END
#SBATCH -p cclake

#SBATCH -o logs/bamMetrics-%j.out

. /etc/profile.d/modules.sh                # Leave this line (enables the module command)
module purge                               # Removes all modules still loaded
module load rhel8/default-ccl              # REQUIRED - loads the basic environment


SAMPLE=$1
REF=$2
source ${SAMPLE}.config

module load $PICARD

BAM_FILE=${SAMPLE}-${REF}.bam
STATS_FILE="${SAMPLE}-${REF}.stats"


${SAMTOOLS} flagstat ${BAM_FILE} > metrics/${SAMPLE}.flagstat
picard CollectWgsMetrics INPUT=${BAM_FILE} O=metrics/${SAMPLE}.wgs.txt REFERENCE_SEQUENCE=${FASTA}/${GENOME}.fasta STOP_AFTER=100000000 VALIDATION_STRINGENCY=LENIENT
picard CollectInsertSizeMetrics INPUT=${BAM_FILE} O=metrics/${SAMPLE}.insert_size.txt H=metrics/${SAMPLE}.insert_size_histogram.pdf VALIDATION_STRINGENCY=LENIENT

echo ${SAMPLE}
echo "..."

grep 'in total' metrics/${SAMPLE}.flagstat
grep '%' metrics/${SAMPLE}.flagstat | grep -v 'singleton'
echo "..."

head -8 metrics/${SAMPLE}.insert_size.txt | tail -2 | cut -f 6,8
echo "..."

head -8 metrics/${SAMPLE}.wgs.txt | tail -2 | cut -f 2,3,15,16,18,20-22
echo "..."


echo -e "BAM\t"`date -r ${BAM_FILE} +"%Y-%m-%d %H:%M:%S"` >> ${STATS_FILE}
echo -e "bamMd5\t"`cat ${SAMPLE}-${REF}.bam.md5` >> ${STATS_FILE}

echo -e "totalReads\t"`grep 'in total' metrics/${SAMPLE}.flagstat | cut -f 1 -d ' '` >> ${STATS_FILE}
echo -e "primaryReads\t"`grep 'primary' metrics/${SAMPLE}.flagstat | grep -v mapped | grep -v duplicates | cut -f 1 -d ' '` >> ${STATS_FILE}
echo -e "mappedReads\t"`grep 'mapped (' metrics/${SAMPLE}.flagstat | grep -v 'primary' | cut -f 1 -d ' '` >> ${STATS_FILE}
echo -e "primaryMapped\t"`grep 'primary mapped' metrics/${SAMPLE}.flagstat | cut -f 1 -d ' '` >> ${STATS_FILE}
echo -e "pairedReads\t"`grep 'properly paired' metrics/${SAMPLE}.flagstat | cut -f 1 -d ' '` >> ${STATS_FILE}

echo -e "meanInsertSize\t"`head -8 metrics/${SAMPLE}.insert_size.txt | tail -1 | cut -f 6` >> ${STATS_FILE}
echo -e "meanReadDepth\t"`head -8 metrics/${SAMPLE}.wgs.txt | tail -1 | cut -f 2` >> ${STATS_FILE}