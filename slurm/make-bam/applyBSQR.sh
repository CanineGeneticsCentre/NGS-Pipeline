#!/usr/bin/env bash

#! RUN : sbatch applyBSQR.sh <SAMPLE>

#! sbatch directives begin here ###############################
#! How many whole nodes should be allocated?
#SBATCH --nodes=1
#! How many (MPI) tasks will there be in total? (<= nodes*32)
#! The skylake/skylake-himem nodes have 32 CPUs (cores) each.
#SBATCH --ntasks=1
#! How much wallclock time will be required?
#SBATCH --time 04:00:00
#! What types of email messages do you wish to receive?
#SBATCH --mail-type=FAIL,INVALID_DEPEND
#! Uncomment this to prevent the job from being requeued (e.g. if
#! interrupted by node failure or system downtime):
##SBATCH --no-requeue
#! For 6GB per CPU, set "-p skylake"; for 12GB per CPU, set "-p skylake-himem":
#SBATCH -p skylake
#SBATCH --mem=5gb

#SBATCH -o logs/applyBQSR-%A_%a.out

. /etc/profile.d/modules.sh                 # Leave this line (enables the module command)
module purge                                # Removes all modules still loaded
module load rhel7/default-peta4             # REQUIRED - loads the basic environment


SAMPLE=$1
source ${SAMPLE}.config
n=$(printf "%04d" $SLURM_ARRAY_TASK_ID)

module load ${GATK}

gatk --java-options "-Djava.io.tmpdir=${HOME}/hpc-work/tmp/ -Xmx2G" ApplyBQSR \
  -R ${FASTA}/${GENOME}.fasta \
  -I ${SAMPLE}.sorted.bam \
  -O base_recal/${SAMPLE}.${n}.bam \
  -L intervals/${n}-scattered.interval_list \
  -bqsr ${SAMPLE}.bsqr.out \
  --static-quantized-quals 10 \
  --static-quantized-quals 20 \
  --static-quantized-quals 30 \
  --add-output-sam-program-record \
  --create-output-bam-md5 \
  --use-original-qualities 
