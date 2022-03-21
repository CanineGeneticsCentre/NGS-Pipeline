#!/usr/bin/env bash

#! RUN : sbatch gvcf2GenomicsDB.sh <SAMPLE_LIST> <REF>

#! sbatch directives begin here ###############################
#! How many whole nodes should be allocated?
#SBATCH --nodes=1
#! How many (MPI) tasks will there be in total? (<= nodes*32)
#! The skylake/skylake-himem nodes have 32 CPUs (cores) each.
#SBATCH --ntasks=8
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


module purge                                # Removes all modules still loaded
module load rhel7/default-peta4             # REQUIRED - loads the basic environment

#module load jdk-8u141-b15-gcc-5.4.0-p4aaopt 
#module load gatk/4.1.0.0
module load gatk-4.2.5.0-gcc-5.4.0-hzdcjga

SAMPLE_LIST=$1
REF=$2
source ${REF}.config


INTERVALS=`head -${SLURM_ARRAY_TASK_ID} ${FASTA}/genomicsDB.intervals | tail -1 | sed s/" "/" -L "/g`
CHR=`echo ${INTERVALS} | cut -f 1 -d' ' | cut -d'_' -f 1`

if [[ ${#CHR} -lt 4 ]] ; then
  CHR="chr"${CHR}
fi


GVCFs=""
for s in `cat ${SAMPLE_LIST}`; do GVCFs+="-V ${s}-${REF}.g.vcf.gz "; done

rm -rf ${GDB}/${GENOME}/${CHR}

gatk --java-options "-Djava.io.tmpdir=${HOME}/hpc-work/tmp/ -Xmx10G" GenomicsDBImport \
    ${GVCFs} \
    --tmp-dir=${HOME}/hpc-work/tmp/ \
    --genomicsdb-workspace-path ${GDB}/${GENOME}/${CHR} \
    -L ${INTERVALS} 
#    --overwrite-existing-genomicsdb-workspace
#    --batch-size 50
#    --reader-threads 8
