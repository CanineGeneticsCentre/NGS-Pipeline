#!/usr/bin/env bash

#! RUN : sbatch annotateVcf.sh <REF>

#! sbatch directives begin here ###############################
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time 02:00:00
#SBATCH --mail-type=FAIL,INVALID_DEPEND,END
##SBATCH --no-requeue
#SBATCH -p cclake

#SBATCH -o logs/snpeff-%j.out

. /etc/profile.d/modules.sh                # Leave this line (enables the module command)
module purge                               # Removes all modules still loaded
module load rhel8/default-ccl              # REQUIRED - loads the basic environment

REF=$1
CHR=$2

source ${REF}.config

module load ${TABIX}

#ls -1v ${CHR}/*.filtered.vcf.gz > ${CHR}.list
#bcftools concat -f ${CHR}.list > ${REF}-${CHR}.vcf

${SNPEFF} -v -csvStats snpEff/${REF}-${CHR}.csv \
    -stats snpEff/${REF}-${CHR}.html \
    ${SNPEFF_DB} \
    ${REF}-${CHR}.filtered.vcf.gz | ${BGZIP} -c > ${REF}-${CHR}.ann.vcf.gz

${TABIX} -p vcf ${REF}-${CHR}.ann.vcf.gz

#if [ -s "${REF}-${CHR}.filtered.vcf.gz" ]; then rm -rf ${REF}-${CHR}.filtered.vcf.gz ${REF}-${CHR}.filtered.vcf.gz.tbi; fi

if [ ! -f samples.list ]
then
  touch samples.list
  for s in `${BCFTOOLS} query -l ${REF}-${CHR}.ann.vcf.gz`; do 
    echo -e "$s\tOmit" >> samples.list
  done
fi
