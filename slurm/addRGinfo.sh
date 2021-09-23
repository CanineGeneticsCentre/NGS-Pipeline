#!/usr/bin/env bash

#! RUN : sbatch sam2bam.sh <SAMPLE>>

#! sbatch directives begin here ###############################
#! How many whole nodes should be allocated?
#SBATCH --nodes=1
#! How many (MPI) tasks will there be in total? (<= nodes*32)
#! The skylake/skylake-himem nodes have 32 CPUs (cores) each.
#SBATCH --ntasks=1
#! How much wallclock time will be required?
#SBATCH --time 01:00:00
#! What types of email messages do you wish to receive?
#SBATCH --mail-type=ALL
#! Uncomment this to prevent the job from being requeued (e.g. if
#! interrupted by node failure or system downtime):
##SBATCH --no-requeue
#! For 6GB per CPU, set "-p skylake"; for 12GB per CPU, set "-p skylake-himem":
#SBATCH -p skylake

#SBATCH -o ../logs/job-%A_%a.out

module purge                                # Removes all modules still loaded
module load rhel7/default-peta4             # REQUIRED - loads the basic environment

#module load samtools/1.10                   # samtools
#module load picard/2.9.2                    # picard

module load jdk-8u141-b15-gcc-5.4.0-p4aaopt 
module load gatk/4.1.0.0


SAMPLE=$1
source ${SAMPLE}.config

#export PICARD_JAVA_MEM_MX='12g'
#export PICARD_JAVA_TMPDIR="${HOME}/rds/hpc-work/tmp"

#picard_latest AddOrReplaceReadGroups INPUT=${SAMPLE}.s_${SLURM_ARRAY_TASK_ID}.aligned.bam OUTPUT=${SAMPLE}.s_${SLURM_ARRAY_TASK_ID}_aligned_sorted_rg.bam RGID=${BARCODE}.${SLURM_ARRAY_TASK_ID} RGLB=${SAMPLE} RGPL=${PLATFORM} RGPU=${BARCODE}.${SLURM_ARRAY_TASK_ID} RGSM=${SAMPLE} SORT_ORDER=coordinate CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT


gatk AddOrReplaceReadGroups --INPUT ${SAMPLE}.s_${SLURM_ARRAY_TASK_ID}.aligned.bam --OUTPUT ${SAMPLE}.s_${SLURM_ARRAY_TASK_ID}_aligned_rg.bam --RGLB ${SAMPLE} --RGPL ${PLATFORM} --RGPU ${BARCODE}.${SLURM_ARRAY_TASK_ID} --RGSM ${SAMPLE} --RGID ${BARCODE}.${SLURM_ARRAY_TASK_ID} --TMP_DIR ${HOME}/hpc-work/tmp/ --VALIDATION_STRINGENCY SILENT

bam_size=$(wc -c < ${SAMPLE}.s_${SLURM_ARRAY_TASK_ID}_aligned_rg.bam)
if [ $bam_size -ge 50000000 ]; then
	rm -rf ${SAMPLE}.s_${SLURM_ARRAY_TASK_ID}.aligned.bam
fi

## NOTES on read groups....
#ID: sample_name.flowcell.lane.barcode
#SM: sample_name
#PL: technology, i.e. ILLUMINA
#PU: flowcell.lane
#LB: sample_name.library_preparation

#ID = Read group identifier = {FLOWCELL_BARCODE}.{LANE}
#PU = Platform Unit = {FLOWCELL_BARCODE}.{LANE}.{library-specific identifier}. This is the most specific definition for a group of reads.
#SM = biological sample name.
#PL = Platform (Valid values: ILLUMINA, SOLID, LS454, HELICOS and PACBIO)
#LB = name of DNA preparation library tube = {SM}.{library-specific identifier}(Important To identify PCR duplicates in MarkDuplicates step. Ignore in PCR free libraries)
