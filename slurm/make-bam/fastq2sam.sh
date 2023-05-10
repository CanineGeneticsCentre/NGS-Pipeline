#!/usr/bin/env bash

#! sbatch directives begin here ###############################
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time 02:00:00
#SBATCH --mail-type=BEGIN,END,FAIL,INVALID_DEPEND
#SBATCH -p skylake

#SBATCH -o logs/fastq2sam_%A-%a.out

. /etc/profile.d/modules.sh                 # Leave this line (enables the module command)
module purge                                # Removes all modules still loaded
module load rhel7/default-peta4             # REQUIRED - loads the basic environment


SAMPLE=$1
LANE=$SLURM_ARRAY_TASK_ID
DIR=lane${LANE}

source ${SAMPLE}.config
mkdir $DIR

module load ${GATK}

FILES=(../*.s_${LANE}.*fq.gz)
BARCODE=$(echo ${FILES[0]} | sed 's/..\///' | cut -d'.' -f 2)	#UDP0021
FLOWCELL=$(echo ${FILES[0]} | sed 's/..\///' | cut -d'.' -f 3)	#HGM3FDSX2

gatk --java-options "-Djava.io.tmpdir=${HOME}/hpc-work/tmp/ -Xmx4G" FastqToSam \
  --FASTQ ${FILES[0]} \
  --FASTQ2 ${FILES[1]} \
  --OUTPUT ${DIR}/${SAMPLE}.L${LANE}.unaligned.bam \
  --READ_GROUP_NAME ${FLOWCELL}.L${LANE} \
  --SAMPLE_NAME ${SAMPLE} \
  --LIBRARY_NAME ${SAMPLE} \
  --PLATFORM_UNIT ${FLOWCELL}.L${LANE}.${BARCODE} \
  --PLATFORM illumina \
  --SEQUENCING_CENTER CRUK-CI 
#  --RUN_DATE 2021-06-11T00:00:00-0400
