#!/usr/bin/env bash

#! RUN : sbatch filterVcf.sh <SAMPLE> REF>

#! sbatch directives begin here ###############################
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --time 00:30:00
#SBATCH --mail-type=FAIL
#SBATCH -p cclake-himem
#SBATCH --mem=12000

#SBATCH -o logs/filterVCF_%j.out

. /etc/profile.d/modules.sh                # Leave this line (enables the module command)
module purge                               # Removes all modules still loaded
module load rhel8/default-ccl              # REQUIRED - loads the basic environment

SAMPLE=$1
REF=$2
source ${SAMPLE}.config

module load ${GATK}

gatk --java-options "-Djava.io.tmpdir=${HOME}/hpc-work/tmp/ -Xmx10G" VariantFiltration \
    -R ${FASTA}/${GENOME}.fasta \
    -V ${SAMPLE}-${REF}.vcf.gz \
    --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || SOR > 3.0 || QUAL < 30" \
    --filter-name "basic" \
    -O ${SAMPLE}-${REF}.filtered.vcf.gz

${BCFTOOLS} view -v indels ${SAMPLE}-${REF}.filtered.vcf.gz | bcftools annotate -x INFO,^FORMAT/GT,FORMAT/AD,FORMAT/DP,FORMAT/GQ,FORMAT/PL -Oz -o ${SAMPLE}-${REF}.indels.vcf.gz
${BCFTOOLS} view -v snps ${SAMPLE}-${REF}.filtered.vcf.gz | bcftools annotate -x INFO,^FORMAT/GT,FORMAT/AD,FORMAT/DP,FORMAT/GQ,FORMAT/PL -Oz -o ${SAMPLE}-${REF}.snps.vcf.gz

echo -e "VCF\t"`date -r ${SAMPLE}-${REF}.filtered.vcf.gz +"%Y-%m-%d %H:%M:%S"` >> ${SAMPLE}-${REF}.stats