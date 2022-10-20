#!/usr/bin/env bash

#! RUN : sbatch haplotypeCaller.sh <SAMPLE>

#! sbatch directives begin here ###############################
#! How many whole nodes should be allocated?
#SBATCH --nodes=1
#! How many (MPI) tasks will there be in total? (<= nodes*32)
#! The skylake/skylake-himem nodes have 32 CPUs (cores) each.
#SBATCH --ntasks=1
#! How much wallclock time will be required?
#SBATCH --time 06:00:00
#! What types of email messages do you wish to receive?
#SBATCH --mail-type=FAIL,INVALID_DEPEND
#! Uncomment this to prevent the job from being requeued (e.g. if
#! interrupted by node failure or system downtime):
##SBATCH --no-requeue
#! For 6GB per CPU, set "-p skylake"; for 12GB per CPU, set "-p skylake-himem":
#SBATCH -p skylake
#SBATCH --mem=5gb

#SBATCH -o ../logs/job-%A_%a.out

. /etc/profile.d/modules.sh                 # Leave this line (enables the module command)
module purge                                # Removes all modules still loaded
module load rhel7/default-peta4             # REQUIRED - loads the basic environment

module load jdk-8u141-b15-gcc-5.4.0-p4aaopt 
module load gatk/4.1.0.0                    # GATK 4.1

SAMPLE=$1
REF=$2
PCR_MODEL=$3
source ${SAMPLE}.config

intervals=`head -${SLURM_ARRAY_TASK_ID} sequence_grouping_with_unmapped.txt | tail -1 | sed s/"\t"/" -L "/g`

gatk --java-options "-Djava.io.tmpdir=${HOME}/hpc-work/tmp/ -Xmx4G" HaplotypeCaller \
  -R ${FASTA}/${GENOME}.fasta \
  -I ${SAMPLE}-${REF}.bam \
  -L ${intervals} \
  -O ${SAMPLE}-${REF}.${SLURM_ARRAY_TASK_ID}.g.vcf \
  -GQB 10 -GQB 20 -GQB 30 -GQB 40 -GQB 50 -GQB 60 -GQB 70 -GQB 80 -GQB 90 \
  -ERC GVCF \
  --pcr-indel-model ${PCR_MODEL}
