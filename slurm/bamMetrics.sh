#!/usr/bin/env bash

#! RUN : sbatch bamMetrics.sh <SAMPLE> <REF>

#! sbatch directives begin here ###############################
#! How many whole nodes should be allocated?
#SBATCH --nodes=1
#! How many (MPI) tasks will there be in total? (<= nodes*32)
#! The skylake/skylake-himem nodes have 32 CPUs (cores) each.
#SBATCH --ntasks=1
#! How much wallclock time will be required?
#SBATCH --time 12:00:00
#! What types of email messages do you wish to receive?
#SBATCH --mail-type=ALL
#! Uncomment this to prevent the job from being requeued (e.g. if
#! interrupted by node failure or system downtime):
##SBATCH --no-requeue
#! For 6GB per CPU, set "-p skylake"; for 12GB per CPU, set "-p skylake-himem":
#SBATCH -p skylake
#SBATCH --mem=5gb

#SBATCH -o ../logs/job-%j.out

module purge                                # Removes all modules still loaded
module load rhel7/default-peta4             # REQUIRED - loads the basic environment

module load picard/2.9.2                    # Picard
module load samtools/1.10                   # Samtools

SAMPLE=$1
REF=$2
source ${SAMPLE}.config

BAM_FILE=${SAMPLE}-${REF}.bam


samtools flagstat ${BAM_FILE} > flagstat.out
picard_latest CollectWgsMetrics INPUT=${BAM_FILE} O=collect_wgs_metrics.txt REFERENCE_SEQUENCE=${FASTA}/${GENOME}.fasta STOP_AFTER=100000000 VALIDATION_STRINGENCY=LENIENT
picard_latest CollectInsertSizeMetrics INPUT=${BAM_FILE} O=insert_size_metrics.txt H=insert_size_histogram.pdf VALIDATION_STRINGENCY=LENIENT

echo ${SAMPLE}
echo "..."

grep 'in total' flagstat.out
grep '%' flagstat.out | grep -v 'singleton'
echo "..."

head -8 collect_wgs_metrics.txt | tail -2 | cut -f 2,3,13-18
echo "..."

head -8 insert_size_metrics.txt | tail -2 | cut -f 6,8
echo "..."
