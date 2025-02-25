#!/usr/bin/env bash

#! RUN : sbatch filterMultiVcf.sh -A ${ACCOUNT} -J ${REF}.filterVCF --array=1,30-31

#! sbatch directives begin here ###############################
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time 00:20:00
#SBATCH --mail-type=FAIL,INVALID_DEPEND,END
##SBATCH --no-requeue
#SBATCH -p cclake

#SBATCH -o logs/filter-%A_%a.out

. /etc/profile.d/modules.sh                # Leave this line (enables the module command)
module purge                               # Removes all modules still loaded
module load rhel8/default-ccl              # REQUIRED - loads the basic environment

REF=$1
CHR=$2

source ${REF}.config
module load ${GATK}

gatk --java-options "-Djava.io.tmpdir=${HOME}/hpc-work/tmp/ -Xmx3G" VariantFiltration \
    -R ${FASTA}/${GENOME}.fasta \
    -V ${CHR}/${REF}-${CHR}-${SLURM_ARRAY_TASK_ID}.vcf.gz \
    --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || SOR > 3.0 || QUAL < 30" \
    --filter-name "basic" \
    -O ${CHR}/${REF}-${CHR}-${SLURM_ARRAY_TASK_ID}.filtered.vcf.gz

#if [ -s "${CHR}/${REF}-${CHR}-${SLURM_ARRAY_TASK_ID}.vcf.gz" ]; then
#  mv ${CHR}/${REF}-${CHR}-${SLURM_ARRAY_TASK_ID}.vcf.gz ${CHR}/${REF}-${CHR}-${SLURM_ARRAY_TASK_ID}.vcf.gz.tbi ${CHR}/done
#  #rm -rf ${CHR}/${REF}-${CHR}-${SLURM_ARRAY_TASK_ID}.vcf.gz ${CHR}/${REF}-${CHR}-${SLURM_ARRAY_TASK_ID}.vcf.gz.tbi
#fi