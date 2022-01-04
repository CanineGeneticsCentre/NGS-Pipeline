#!/usr/bin/env bash

#! RUN : sbatch mergeGvcf.sh <SAMPLE> <INTERVALS> <REF>

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

module load jdk-8u141-b15-gcc-5.4.0-p4aaopt 
module load gatk/4.1.0.0                    # GATK 4.1

SAMPLE=$1
INTERVALS=$2
REF=$3
source ${SAMPLE}.config

#ls -1 ${SAMPLE}*.recal_data.csv > bsqr_reports.txt
GVFS=""
for i in `seq 1 ${INTERVALS}`; do GVFS+="-I ${SAMPLE}.$i.gvcf "; done

gatk --java-options  "-Djava.io.tmpdir=${HOME}/hpc-work/tmp/ -Xmx2G" MergeVcfs ${GVFS} -O ${SAMPLE}-${REF}.gvcf

#gvcf_size=$(wc -c < ${SAMPLE}-${REF}.gvcf)
#if [ $gvcf_size -ge 50000000 ];then
#  for i in `seq 1 ${INTERVALS}`; do 
#	  rm -rf ${SAMPLE}.$i.gvcf
#  done
#fi
