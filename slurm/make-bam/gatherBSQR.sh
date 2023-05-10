#!/usr/bin/env bash

#! RUN : sbatch gatherBSQR.sh <SAMPLE>

#! sbatch directives begin here ###############################
#! How many whole nodes should be allocated?
#SBATCH --nodes=1
#! How many (MPI) tasks will there be in total? (<= nodes*32)
#! The skylake/skylake-himem nodes have 32 CPUs (cores) each.
#SBATCH --ntasks=1
#! How much wallclock time will be required?
#SBATCH --time 01:00:00
#! What types of email messages do you wish to receive?
#SBATCH --mail-type=FAIL,INVALID_DEPEND
#! Uncomment this to prevent the job from being requeued (e.g. if
#! interrupted by node failure or system downtime):
##SBATCH --no-requeue
#! For 6GB per CPU, set "-p skylake"; for 12GB per CPU, set "-p skylake-himem":
#SBATCH -p skylake
#SBATCH --mem=5gb

#SBATCH -o logs/gatherBQSR-%j.out

. /etc/profile.d/modules.sh                 # Leave this line (enables the module command)
module purge                                # Removes all modules still loaded
module load rhel7/default-peta4             # REQUIRED - loads the basic environment


SAMPLE=$1
source ${SAMPLE}.config

module load ${GATK}

REPORTS=""
for f in `ls base_recal/${SAMPLE}.*.bsqr.txt`; do REPORTS+="-I $f "; done


gatk --java-options  "-Djava.io.tmpdir=${HOME}/hpc-work/tmp/ -Xmx4G" GatherBQSRReports ${REPORTS} -O ${SAMPLE}.bsqr.out

if [ -s "${SAMPLE}.bsqr.out" ]; then
	rm -rf base_recal/${SAMPLE}.*.bsqr.txt
fi
