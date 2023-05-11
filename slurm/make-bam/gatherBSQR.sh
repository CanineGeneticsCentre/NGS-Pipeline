#!/usr/bin/env bash

#! RUN : sbatch gatherBSQR.sh <SAMPLE>

#! sbatch directives begin here ###############################
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time 00:30:00
#SBATCH --mail-type=FAIL,INVALID_DEPEND
##SBATCH --no-requeue
#SBATCH -p cclake

#SBATCH -o logs/gatherBQSR-%j.out

. /etc/profile.d/modules.sh                 # Leave this line (enables the module command)
module purge
module load rhel7/default-ccl


SAMPLE=$1
source ${SAMPLE}.config

module load ${GATK}

REPORTS=""
for f in `ls base_recal/${SAMPLE}.*.bsqr.txt`; do REPORTS+="-I $f "; done


gatk --java-options  "-Djava.io.tmpdir=${HOME}/hpc-work/tmp/ -Xmx2G" GatherBQSRReports ${REPORTS} -O ${SAMPLE}.bsqr.out

if [ -s "${SAMPLE}.bsqr.out" ]; then
	rm -rf base_recal/${SAMPLE}.*.bsqr.txt
fi
