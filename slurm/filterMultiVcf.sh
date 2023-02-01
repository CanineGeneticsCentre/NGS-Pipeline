#!/usr/bin/env bash

#! RUN : sbatch filterMultiVcf.sh

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
#SBATCH -p skylake-himem
#SBATCH --mem=10gb

#SBATCH -o logs/job-%A_%a.out

. /etc/profile.d/modules.sh                 # Leave this line (enables the module command)
module purge                                # Removes all modules still loaded
module load rhel7/default-peta4             # REQUIRED - loads the basic environment

module load gatk-4.2.5.0-gcc-5.4.0-hzdcjga

source ${REF}.config

INTERVALS=`head -${SLURM_ARRAY_TASK_ID} ${FASTA}/${REF}-genomicsDB.intervals | tail -1 | sed s/" "/" -L "/g`
CHR=`echo ${INTERVALS} | cut -f 1 -d' ' | cut -d'_' -f 1 | cut -f 1 -d':'`

if [[ ${#CHR} -lt 4 ]] ; then
  CHR="chr"${CHR}
fi

gatk --java-options "-Djava.io.tmpdir=${HOME}/hpc-work/tmp/ -Xmx10G" VariantFiltration \
    -R ${FASTA}/${GENOME}.fasta \
    -V ${CHR}/${REF}-${CHR}-${SLURM_ARRAY_TASK_ID}.vcf.gz \
    --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || SOR > 3.0 || QUAL < 30" \
    --filter-name "basic" \
    -O ${CHR}/${REF}-${CHR}-${SLURM_ARRAY_TASK_ID}.filtered.vcf.gz

rm -rf ${CHR}/${REF}-${CHR}-${SLURM_ARRAY_TASK_ID}.vcf.gz ${CHR}/${REF}-${CHR}-${SLURM_ARRAY_TASK_ID}.vcf.gz.tbi