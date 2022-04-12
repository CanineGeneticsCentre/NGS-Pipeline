#!/usr/bin/env bash

#! RUN : sbatch filterVcf.sh <REF>

#! sbatch directives begin here ###############################
#! How many whole nodes should be allocated?
#SBATCH --nodes=1
#! How many (MPI) tasks will there be in total? (<= nodes*32)
#! The skylake/skylake-himem nodes have 32 CPUs (cores) each.
#SBATCH --ntasks=1
#! How much wallclock time will be required?
#SBATCH --time 00:30:00
#! What types of email messages do you wish to receive?
#SBATCH --mail-type=FAIL,INVALID_DEPEND,END
#! Uncomment this to prevent the job from being requeued (e.g. if
#! interrupted by node failure or system downtime):
##SBATCH --no-requeue
#! For 6GB per CPU, set "-p skylake"; for 12GB per CPU, set "-p skylake-himem":
#SBATCH -p skylake-himem
#SBATCH --mem=10gb

#SBATCH -o logs/job-%j.out


module purge                                # Removes all modules still loaded
module load rhel7/default-peta4             # REQUIRED - loads the basic environment

module load gatk-4.2.5.0-gcc-5.4.0-hzdcjga

REF=$1
CHR=$2
ID=$3
source ${REF}.config

gatk --java-options "-Djava.io.tmpdir=${HOME}/hpc-work/tmp/ -Xmx10G" VariantFiltration \
    -R ${FASTA}/${GENOME}.fasta \
    -V ${CHR}/${REF}-${CHR}-${ID}.vcf.gz \
    --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || SOR > 3.0 || QUAL < 30" \
    --filter-name "basic" \
    -O ${CHR}/${REF}-${CHR}-${ID}.filtered.vcf.gz

#rm -rf ${CHR}/${REF}-${CHR}-${ID}.vcf.gz ${CHR}/${REF}-${CHR}-${ID}.vcf.gz.tbi