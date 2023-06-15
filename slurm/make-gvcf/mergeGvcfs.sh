#!/usr/bin/env bash

#! RUN : sbatch mergeGvcfs.sh <SAMPLE>

#! sbatch directives begin here ###############################
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time 06:00:00
#SBATCH --mail-type=FAIL,END
##SBATCH --no-requeue
#SBATCH -p cclake

#SBATCH -o logs/gvcfMerge-%j.out

. /etc/profile.d/modules.sh                 # Leave this line (enables the module command)
module purge
module load rhel7/default-ccl


SAMPLE=$1
source ${SAMPLE}.config

module load ${GATK}

GVCFS=""
for i in `seq 0 $(($INTERVALS-1))`; do n=$(printf "%04d" $i); GVCFS+="--INPUT gvcf/${SAMPLE}-${REF}.$n.g.vcf "; done


gatk --java-options "-Djava.io.tmpdir=${HOME}/hpc-work/tmp/ -Xmx3G" MergeVcfs \
  ${GVCFS} \
  -O ${SAMPLE}-${REF}.g.vcf.gz

gvcf_size=$(wc -c < ${SAMPLE}-${REF}.g.vcf.gz)
if [ $gvcf_size -ge 1000000 ]; then
  rm -rf gvcf
fi
