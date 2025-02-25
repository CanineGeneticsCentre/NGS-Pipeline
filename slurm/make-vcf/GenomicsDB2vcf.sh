#!/usr/bin/env bash

#! RUN : sbatch gvcf2GenomicsDB.sh

#! sbatch directives begin here ###############################
#SBATCH --ntasks=1
#SBATCH --time 36:00:00
#SBATCH --mail-type=ALL
##SBATCH --no-requeue
#SBATCH -p cclake-himem
#SBATCH --cpus-per-task 4   ## hopefully gives me 4 cores on the same node, therefore 4 x 6840MB of RAM... (maybe!)

#SBATCH -o logs/gdbExport-%A_%a.out

. /etc/profile.d/modules.sh                # Leave this line (enables the module command)
module purge                               # Removes all modules still loaded
module load rhel8/default-ccl              # REQUIRED - loads the basic environment

REF=$1
CHR=$2

source ${REF}.config
module load ${GATK}

# one job array per chr; chunk chr into bits [chrX.intervals] for each job array part.
LOCATION=`head -${SLURM_ARRAY_TASK_ID} ${CHR}.intervals | tail -1`

cd ${CHR}; mkdir logs


gatk --java-options "-Djava.io.tmpdir=${HOME}/hpc-work/tmp/ -Xmx25G" GenotypeGVCFs \
    -R ${FASTA}/${GENOME}.fasta \
    --tmp-dir ${HOME}/hpc-work/tmp/ \
    -V gendb://${GDB}/${CHR} \
    -L ${LOCATION} \
    -O ${REF}-${CHR}-${SLURM_ARRAY_TASK_ID}.vcf.gz

exit;

#INTERVALS=`head -${SLURM_ARRAY_TASK_ID} ${FASTA}/${REF}-genomicsDB.intervals | tail -1 | sed s/" "/" -L "/g`
#CHR=`echo ${INTERVALS} | cut -f 1 -d' ' | cut -d'_' -f 1 | cut -f 1 -d':'`

#gatk --java-options "-Djava.io.tmpdir=${HOME}/hpc-work/tmp/ -Xmx50G" GenotypeGVCFs \
#    -R ${FASTA}/${GENOME}.fasta \
#    --tmp-dir ${HOME}/hpc-work/tmp/ \
#    -V gendb://${GDB}/${GENOME}/${CHR} \
#    -O ${CHR}/${REF}-${CHR}.vcf.gz
