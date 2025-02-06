#!/usr/bin/env bash

#! sbatch directives begin here ###############################
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time 02:00:00
#SBATCH --mail-type=BEGIN,FAIL
#SBATCH -p cclake-himem

#SBATCH -o logs/fastq2sam_%A-%a.out

. /etc/profile.d/modules.sh                # Leave this line (enables the module command)
module purge                               # Removes all modules still loaded
module load rhel8/default-ccl              # REQUIRED - loads the basic environment


SAMPLE=$1
DIR=rg${SLURM_ARRAY_TASK_ID}

source ${SAMPLE}.config
mkdir $DIR

module load ${GATK}

FQ=$(head -${SLURM_ARRAY_TASK_ID} ../fastq.files | tail -1 | sed s/.r_1.fq.gz//)
FILES=(../${FQ}.*fq.gz)
FASTQ=$(zcat ${FILES[0]} | head -1)
LIBRARY=$(echo ${FILES[0]} | sed 's/..\///' | cut -d'.' -f 1) #SLX-25208
BARCODE=$(echo ${FILES[0]} | sed 's/..\///' | cut -d'.' -f 2)	#UDP0021
#FLOWCELL=$(echo ${FILES[0]} | sed 's/..\///' | cut -d'.' -f 3)	#HGM3FDSX2
FLOWCELL=$(echo $FASTQ | cut -f 3 -d':')  #22FFK7LT4
LANE=$(echo $FASTQ | cut -f 4 -d':')

gatk --java-options "-Djava.io.tmpdir=${HOME}/hpc-work/tmp/ -Xmx4G" FastqToSam \
  --FASTQ ${FILES[0]} \
  --FASTQ2 ${FILES[1]} \
  --OUTPUT ${DIR}/${SAMPLE}.RG${SLURM_ARRAY_TASK_ID}.unaligned.bam \
  --READ_GROUP_NAME ${FLOWCELL}.L${LANE} \
  --SAMPLE_NAME ${SAMPLE} \
  --LIBRARY_NAME ${LIBRARY} \
  --PLATFORM_UNIT ${FLOWCELL}.L${LANE}.${BARCODE} \
  --PLATFORM illumina \
  --SEQUENCING_CENTER CRUK-CI 
#  --RUN_DATE 2021-06-11T00:00:00-0400
