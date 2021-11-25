#!/usr/bin/env bash

#! RUN : sbatch mergeBam.sh <SAMPLE> <INTERVALS>

#! sbatch directives begin here ###############################
#! How many whole nodes should be allocated?
#SBATCH --nodes=1
#! How many (MPI) tasks will there be in total? (<= nodes*32)
#! The skylake/skylake-himem nodes have 32 CPUs (cores) each.
#SBATCH --ntasks=1
#! How much wallclock time will be required?
#SBATCH --time 04:00:00
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

module load jdk-8u141-b15-gcc-5.4.0-p4aaopt 
module load gatk/4.1.0.0                    # GATK 4.1

SAMPLE=$1
INTERVALS=$2
source ${SAMPLE}.config

#ls -1 ${SAMPLE}*.recal_data.csv > bsqr_reports.txt
BAMS=""
for i in `seq 1 ${INTERVALS}`; do BAMS+="-I ${SAMPLE}.$i.bam "; done

gatk --java-options  "-Djava.io.tmpdir=${HOME}/hpc-work/tmp/ -Xmx2G" GatherBamFiles ${BAMS} -O /dev/stdout | gatk --java-options "-Djava.io.tmpdir=${HOME}/hpc-work/tmp/ -Xmx2G" SortSam --INPUT /dev/stdin  --OUTPUT ${SAMPLE}.bam --SORT_ORDER "coordinate" --TMP_DIR ${HOME}/hpc-work/tmp/ --CREATE_INDEX true --CREATE_MD5_FILE true

#gatk --java-options  "-Djava.io.tmpdir=${HOME}/hpc-work/tmp/ -Xmx2G" GatherBamFiles ${BAMS} -O ${SAMPLE}.merged.bam


bam_size=$(wc -c < ${SAMPLE}.bam)
if [ $bam_size -ge 50000000 ];then
  for i in `seq 1 ${INTERVALS}`; do 
	  rm -rf ${SAMPLE}.$i.bam
	  rm -rf ${SAMPLE}.$i.bai
	  rm -rf ${SAMPLE}.$i.bam.md5
  done
  rm -rf ${SAMPLE}.sorted.bam ${SAMPLE}.sorted.bai ${SAMPLE}.sorted.bam.md5
fi
