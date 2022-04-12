#!/usr/bin/env bash

#! RUN : sbatch filterVcf.sh <REF>

#! sbatch directives begin here ###############################
#! How many whole nodes should be allocated?
#SBATCH --nodes=1
#! How many (MPI) tasks will there be in total? (<= nodes*32)
#! The skylake/skylake-himem nodes have 32 CPUs (cores) each.
#SBATCH --ntasks=1
#! How much wallclock time will be required?
#SBATCH --time 12:00:00
#! What types of email messages do you wish to receive?
#SBATCH --mail-type=FAIL,INVALID_DEPEND,END
#! Uncomment this to prevent the job from being requeued (e.g. if
#! interrupted by node failure or system downtime):
##SBATCH --no-requeue
#! For 6GB per CPU, set "-p skylake"; for 12GB per CPU, set "-p skylake-himem":
#SBATCH -p skylake
#SBATCH --mem=5gb

#SBATCH -o logs/job-%j.out


module purge                                # Removes all modules still loaded
module load rhel7/default-peta4             # REQUIRED - loads the basic environment

module load bcftools-1.9-gcc-5.4.0-b2hdt5n
module load tabix-2013-12-16-gcc-5.4.0-xn3xiv7

CHR=$1
source ${REF}.config

if [[ ${#CHR} -lt 4 ]] ; then
  CHR="chr"${CHR}
fi

ls -1 ${CHR}/*.filtered.vcf.gz > ${CHR}.list

bcftools concat -f ${CHR}.list > ${REF}-${CHR}.vcf

/rds/project/rds-Qr3fy2NTCy0/Software/local/snpEff/scripts/snpEff \
    -v -csvStats snpEff/${REF}-${CHR}.csv \
    -stats snpEff/${REF}-${CHR}.html \
    ${SNPEFF} \
    ${REF}-${CHR}.vcf | bgzip -c > ${REF}-${CHR}.ann.vcf.gz

rm -rf ${REF}-${CHR}.vcf
tabix -p vcf ${REF}-${CHR}.ann.vcf.gz