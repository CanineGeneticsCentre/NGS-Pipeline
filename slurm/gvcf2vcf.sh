#!/usr/bin/env bash

#! RUN : sbatch gvcf2vcf.sh <SAMPLE> <REF>

#! sbatch directives begin here ###############################
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time 04:00:00
#SBATCH --mail-type=FAIL
#SBATCH -p cclake
#SBATCH --mem=5gb

#SBATCH -o logs/makeVCF_%j.out

. /etc/profile.d/modules.sh                # Leave this line (enables the module command)
module purge                               # Removes all modules still loaded
module load rhel8/default-ccl              # REQUIRED - loads the basic environment


SAMPLE=$1
REF=$2

source ${SAMPLE}.config

module load ${GATK}

gatk --java-options "-Djava.io.tmpdir=${HOME}/hpc-work/tmp/ -Xmx4G" GenotypeGVCFs \
    -R ${FASTA}/${GENOME}.fasta \
    -V ${SAMPLE}-${REF}.g.vcf.gz \
    -O ${SAMPLE}-${REF}.vcf.gz
