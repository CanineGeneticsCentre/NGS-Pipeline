#!/usr/bin/env bash

#! RUN : sbatch filterMultiVcf.sh -A ${ACCOUNT} -J ${REF}.filterVCF --array=1,30-31

#! sbatch directives begin here ###############################
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time 01:00:00
#SBATCH --mail-type=FAIL,INVALID_DEPEND,END
##SBATCH --no-requeue
#SBATCH -p cclake

#SBATCH -o logs/concat-%j.out

. /etc/profile.d/modules.sh                # Leave this line (enables the module command)
module purge                               # Removes all modules still loaded
module load rhel8/default-ccl              # REQUIRED - loads the basic environment

REF=$1
CHR=$2

source ${REF}.config
#module load ${BCFTOOLS}
module load ${TABIX}

echo bcftools concat -f ${CHR}/files.list -Oz -o ${REF}-${CHR}.filtered.vcf.gz
bcftools concat -f ${CHR}/files.list -Oz -o ${REF}-${CHR}.filtered.vcf.gz
tabix -p vcf ${REF}-${CHR}.filtered.vcf.gz

#if [ -s "${REF}-${CHR}.filtered.vcf.gz" ]; then rm -rf ${CHR}; fi
