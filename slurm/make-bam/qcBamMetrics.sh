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
#SBATCH --mail-type=FAIL,INVALID_DEPEND,END
#! Uncomment this to prevent the job from being requeued (e.g. if
#! interrupted by node failure or system downtime):
##SBATCH --no-requeue
#! For 6GB per CPU, set "-p skylake"; for 12GB per CPU, set "-p skylake-himem":
#SBATCH -p skylake
#SBATCH --mem=5gb

#SBATCH -o logs/bamMetrics-%j.out

. /etc/profile.d/modules.sh                 # Leave this line (enables the module command)
module purge                                # Removes all modules still loaded
module load rhel7/default-peta4             # REQUIRED - loads the basic environment


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
