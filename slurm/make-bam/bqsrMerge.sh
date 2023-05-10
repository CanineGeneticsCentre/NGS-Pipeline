#!/usr/bin/env bash

#! RUN : sbatch mergeBam.sh <SAMPLE> <INTERVALS>

#! sbatch directives begin here ###############################
#! How many whole nodes should be allocated?
#SBATCH --nodes=1
#! How many (MPI) tasks will there be in total? (<= nodes*32)
#! The skylake/skylake-himem nodes have 32 CPUs (cores) each.
#SBATCH --ntasks=1
#! How much wallclock time will be required?
#SBATCH --time 04:00:00
#! What types of email messages do you wish to receive?
#SBATCH --mail-type=FAIL,INVALID_DEPEND,END
#! Uncomment this to prevent the job from being requeued (e.g. if
#! interrupted by node failure or system downtime):
##SBATCH --no-requeue
#! For 6GB per CPU, set "-p skylake"; for 12GB per CPU, set "-p skylake-himem":
#SBATCH -p skylake
#SBATCH --mem=5gb

#SBATCH -o logs/bqsrMerge-%j.out

. /etc/profile.d/modules.sh                 # Leave this line (enables the module command)
module purge                                # Removes all modules still loaded
module load rhel7/default-peta4             # REQUIRED - loads the basic environment


SAMPLE=$1
INTERVALS=$2
REF=$3
source ${SAMPLE}.config

module load ${GATK}

BAMS=""
for i in `seq 0 $(($INTERVALS-1))`; do n=$(printf "%04d" $i); BAMS+="-I base_recal/${SAMPLE}.$n.bam "; done

gatk --java-options  "-Djava.io.tmpdir=${HOME}/hpc-work/tmp/ -Xmx4G" GatherBamFiles ${BAMS} -O /dev/stdout | \
gatk --java-options "-Djava.io.tmpdir=${HOME}/hpc-work/tmp/ -Xmx4G" SortSam \
  --INPUT /dev/stdin  \
  --OUTPUT ${SAMPLE}-${REF}.bam \
  --SORT_ORDER "coordinate" \
  --TMP_DIR ${HOME}/hpc-work/tmp/ \
  --CREATE_INDEX true \
  --CREATE_MD5_FILE true


if [ ! -f ${SAMPLE}-${REF}.bam ]; then
  echo "ERROR - ${SAMPLE}-${REF}.bam does not exist";
	exit 1;
fi

input_size=$(stat -c%s "${SAMPLE}.sorted.bam")
output_size=$(stat -c%s "${SAMPLE}-${REF}.bam")
if [ output_size > input_size ]; then
  rm -rf base_recal;
  mv ${SAMPLE}.sorted.bam ${SAMPLE}.sorted.bai ${SAMPLE}.sorted.bam.md5 ${SAMPLE}.bsqr.out tmp_files/
fi
