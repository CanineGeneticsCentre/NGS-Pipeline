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

REF=$1
ID=$2
source ${REF}.config

echo /rds/project/rds-Qr3fy2NTCy0/Software/local/snpEff/scripts/snpEff -v ${SNPEFF} ${CHR}/${REF}-${CHR}-${ID}.filtered.vcf.gz > ${CHR}/${REF}-${CHR}-${ID}.final.vcf.gz
/rds/project/rds-Qr3fy2NTCy0/Software/local/snpEff/scripts/snpEff -v ${SNPEFF} ${CHR}/${REF}-${CHR}-${ID}.filtered.vcf.gz > ${CHR}/${REF}-${CHR}-${ID}.final.vcf.gz

#rm -rf ${CHR}/${REF}-${CHR}-${ID}.filtered.vcf.gz ${CHR}/${REF}-${CHR}-${ID}.filtered.vcf.gz.tbi