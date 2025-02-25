#!/usr/bin/env bash

#! sbatch directives begin here ###############################
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time 01:00:00
#SBATCH --mail-type=FAIL
#SBATCH -p cclake

#SBATCH -o logs/count-reads_%j.out

. /etc/profile.d/modules.sh                # Leave this line (enables the module command)
module purge                               # Removes all modules still loaded
module load rhel8/default-ccl              # REQUIRED - loads the basic environment

SAMPLE=$1
REF=$2

for f in `ls ../*R1*.fq.gz ../*.r_1.fq.gz ../*_1.*.fq.gz`; do 
  ((read_count = read_count + $(zcat $f | wc -l))); 
done

((read_count = read_count * 2))
((read_count = read_count / 4))

echo $read_count > metrics/${SAMPLE}-reads.count;
echo -e "rawReads\t"$read_count >> ${SAMPLE}-${REF}.stats