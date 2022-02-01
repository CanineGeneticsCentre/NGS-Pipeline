#!/usr/bin/env bash

#! RUN : sbatch sortSam.sh <SAMPLE>

#! sbatch directives begin here ###############################
#! How many whole nodes should be allocated?
#SBATCH --nodes=1
#! How many (MPI) tasks will there be in total? (<= nodes*32)
#! The skylake/skylake-himem nodes have 32 CPUs (cores) each.
#SBATCH --ntasks=1
#! How much wallclock time will be required?
#SBATCH --time 12:00:00
#! What types of email messages do you wish to receive?
#SBATCH --mail-type=FAIL,INVALID_DEPEND
#! Uncomment this to prevent the job from being requeued (e.g. if
#! interrupted by node failure or system downtime):
##SBATCH --no-requeue
#! For 6GB per CPU, set "-p skylake"; for 12GB per CPU, set "-p skylake-himem":
#SBATCH -p skylake-himem
#SBATCH --mem=10gb

#SBATCH -o ../logs/job-%j.out

module purge                                # Removes all modules still loaded
module load rhel7/default-peta4             # REQUIRED - loads the basic environment

module load jdk-8u141-b15-gcc-5.4.0-p4aaopt 
module load gatk/4.1.0.0                    # GATK 4.1

SAMPLE=$1
source ${SAMPLE}.config


gatk --java-options "-Djava.io.tmpdir=${HOME}/hpc-work/tmp/ -Xmx2G" SortSam  --INPUT ${SAMPLE}.aligned.unsorted.dedup.bam --OUTPUT /dev/stdout --SORT_ORDER "coordinate" --TMP_DIR ${HOME}/hpc-work/tmp/ | gatk --java-options "-Djava.io.tmpdir=${HOME}/hpc-work/tmp/ -Xmx2G" SetNmMdAndUqTags --INPUT /dev/stdin --OUTPUT ${SAMPLE}.sorted.bam --CREATE_INDEX true --CREATE_MD5_FILE true --REFERENCE_SEQUENCE ${FASTA}/${GENOME}.fasta --TMP_DIR ${HOME}/hpc-work/tmp/

bam_size=$(wc -c < ${SAMPLE}.sorted.bam)
if [ $bam_size -ge 50000000 ]; then
	rm -rf ${SAMPLE}.aligned.unsorted.dedup.bam
fi