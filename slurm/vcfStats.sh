#!/bin/bash -l

#SBATCH --export=ALL
#SBATCH --mail-user=ellen.schofield@aht.org.uk
#SBATCH --mail-type=ALL
#SBACTH --mem=1
##SBATCH -e outputs/job-%j.error
#SBATCH -o outputs/job-%j.output
#!/usr/bin/env bash

#! RUN : sbatch vcfStats.sh <LIST>

#! sbatch directives begin here ###############################
#! How many whole nodes should be allocated?
#SBATCH --nodes=1
#! How many (MPI) tasks will there be in total? (<= nodes*32)
#! The skylake/skylake-himem nodes have 32 CPUs (cores) each.
#SBATCH --ntasks=1
#! How much wallclock time will be required?
#SBATCH --time 00:01:00
#! What types of email messages do you wish to receive?
#SBATCH --mail-type=FAIL,END
#! Uncomment this to prevent the job from being requeued (e.g. if
#! interrupted by node failure or system downtime):
##SBATCH --no-requeue
#! For 6GB per CPU, set "-p skylake"; for 12GB per CPU, set "-p skylake-himem":
#SBATCH -p skylake

#SBATCH -o ../logs/job-%j.out

echo
echo ${SAMPLE}
echo "..."

ll ${SAMPLE}-${REF}.vcf.gz
zcat ${SAMPLE}-${REF}.vcf.gz | grep -v '#' | wc -l
echo "..."

ll ${SAMPLE}-${REF}.filtered.vcf.gz
zcat ${SAMPLE}-${REF}.filtered.vcf.gz | grep -v '#' | grep -v 'basic' | wc -l
echo
