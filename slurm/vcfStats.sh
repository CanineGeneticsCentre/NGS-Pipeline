#!/usr/bin/env bash

#! RUN : sbatch vcfStats.sh <LIST>

#! sbatch directives begin here ###############################
#! How many whole nodes should be allocated?
#SBATCH --nodes=1
#! How many (MPI) tasks will there be in total? (<= nodes*32)
#! The skylake/skylake-himem nodes have 32 CPUs (cores) each.
#SBATCH --ntasks=1
#! How much wallclock time will be required?
#SBATCH --time 00:10:00
#! What types of email messages do you wish to receive?
#SBATCH --mail-type=FAIL,END
#! Uncomment this to prevent the job from being requeued (e.g. if
#! interrupted by node failure or system downtime):
##SBATCH --no-requeue
#! For 6GB per CPU, set "-p skylake"; for 12GB per CPU, set "-p skylake-himem":
#SBATCH -p skylake

#SBATCH -o logs/job-%j.out

echo
echo ${SAMPLE} - ${REF}
echo "..."

ls -lh ${SAMPLE}-${REF}.vcf.gz
zcat ${SAMPLE}-${REF}.vcf.gz | grep -v '#' | wc -l
echo "..."

ls -lh ${SAMPLE}-${REF}.filtered.vcf.gz
zcat ${SAMPLE}-${REF}.filtered.vcf.gz | grep -v '#' | grep -v 'basic' | wc -l
echo "..."

echo SNPS
zcat ${SAMPLE}-${REF}.snps.vcf.gz | grep -v '#' | grep -v 'basic' | wc -l
echo Indels
zcat ${SAMPLE}-${REF}.indels.vcf.gz | grep -v '#' | grep -v 'basic' | wc -l
echo


mkdir done
mv ${SAMPLE}-${REF}.vcf.gz* ${SAMPLE}-${REF}.filtered.vcf.gz* done/