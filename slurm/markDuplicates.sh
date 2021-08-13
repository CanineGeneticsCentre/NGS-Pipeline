#!/usr/bin/env bash

#! RUN : sbatch markDuplicates.sh <SAMPLE> <LANES>

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

#SBATCH -o ../logs/job-%j.out

module purge                                # Removes all modules still loaded
module load rhel7/default-peta4             # REQUIRED - loads the basic environment

module load jdk-8u141-b15-gcc-5.4.0-p4aaopt 
module load gatk/4.1.0.0

SAMPLE=$1
LANES=$2
source ../${SAMPLE}.config

INPUT=''
for i in `seq 1 ${LANES}`; do 
  INPUT+=" --INPUT ${SAMPLE}.s_${i}_aligned_rg.bam"
done
echo $INPUT


gatk --java-options "-Xms6G" MarkDuplicates ${INPUT} --OUTPUT ${SAMPLE}.aligned.unsorted.dedup.bam --METRICS_FILE ${SAMPLE}.dedup.metrics --VALIDATION_STRINGENCY=SILENT --ASSUME_SORT_ORDER "queryname"


#  ~{gatk_path} --java-options "-Dsamjdk.compression_level=~{compression_level} -Xms~{command_mem_gb}G" \
#      MarkDuplicates \
#      --INPUT ~{sep=' --INPUT ' input_bams} \
#      --OUTPUT ~{output_bam_basename}.bam \
#      --METRICS_FILE ~{metrics_filename} \
#      --VALIDATION_STRINGENCY SILENT \
#      --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 \
#      --ASSUME_SORT_ORDER "queryname" \
#      --CREATE_MD5_FILE true