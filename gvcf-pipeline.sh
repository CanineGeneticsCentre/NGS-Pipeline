#!/usr/bin/env bash

#! RUN : bash gvcf-pipeline.sh <SAMPLE> <REF>
#! Eg. : bash gvcf-pipeline.sh CS_35365 cf4

SAMPLE=$1
REF=$2
SCRIPTS=`dirname $0`
CFG="${SCRIPTS}/ngs-pipeline-${REF}.config"

[[ -z "$SAMPLE" ]] && { echo "ERROR: No SAMPLE provided for this run"; exit 1; }
[[ -z "$REF" ]] && { echo "ERROR: No REFERENCE provided for this run"; exit 1; }

date

cd $SAMPLE
cp $CFG ${SAMPLE}-${REF}.config; source ${SAMPLE}-${REF}.config

mkdir -p $GENOME; 
mv ${SAMPLE}-${REF}.config $GENOME/${SAMPLE}.config
cd $GENOME

mkdir metrics logs intervals gvcf

#PCR_FREE=false;
PCR_MODEL='CONSERVATIVE';

echo
printf "Is the data PCR free?\n"
printf "\t1. No [default]\n"
printf "\t2. Yes\n"
# Assign input value into a variable
read answer

if [[ -n $answer && $answer == "2" ]]; then
    #PCR_FREE=true;
    PCR_MODEL='NONE';
fi

# Copy BAM/BAI files from RCS... will only copy if files no present or rcs version is newer
rsync --progress -av ${WGS}/${SAMPLE}/${SAMPLE}-${REF}.ba* ./

# generate seqeunce groups for future scatter/gather steps.
#perl ${SCRIPTS}/perl/createSeqGroups.pl ${DICT}
#INTERVALS=`wc -l sequence_grouping.txt | awk '{print $1}'`

module load ${GATK}
gatk SplitIntervals -R ${FASTA}/${GENOME}.fasta -L ${INTERVAL_LIST} --scatter-count ${INTERVALS} -O intervals --subdivision-mode BALANCING_WITHOUT_INTERVAL_SUBDIVISION

#BALANCING_WITHOUT_INTERVAL_SUBDIVISION means that --scatter-count is the upper limit for count of intervals
#Need to reassign INTERVALS to the number actually generated
INTERVALS=`ls intervals | wc -l`
sed -i s/INTERVALS=50/INTERVALS=${INTERVALS}/ ${SAMPLE}.config

# Create gvcf files with HaplotypeCaller
jid1=$(sbatch -A ${ACCOUNT} -J ${SAMPLE}.HC --array=0-$(($INTERVALS-1)) ${SCRIPTS}/slurm/make-gvcf/haplotypeCaller.sh ${SAMPLE} ${REF} ${PCR_MODEL})

# Merge gVCF files into single gVCF
jid2=$(sbatch -A ${ACCOUNT} -J ${SAMPLE}.GVCF --dependency=afterok:${jid1##* } ${SCRIPTS}/slurm/make-gvcf/mergeGvcfs.sh ${SAMPLE})

echo $jid2

###########################################################
#                     VCF Files/Stats                     #
###########################################################

# Convert gVCF file to VCF
jid3=$(sbatch -A ${ACCOUNT} -J ${SAMPLE}.VCF --dependency=afterok:${jid2##* } ${SCRIPTS}/slurm/gvcf2vcf.sh ${SAMPLE} ${REF})

# Filter VCF and split into snps/indels
jid4=$(sbatch -A ${ACCOUNT} -J ${SAMPLE}.filterVCF --dependency=afterok:${jid3##* } ${SCRIPTS}/slurm/filterVcf.sh ${SAMPLE} ${REF})

# Output VCF stats
tmp=$(sbatch -A ${ACCOUNT} -J ${SAMPLE}.VCFstats --dependency=afterok:${jid4##* } ${SCRIPTS}/slurm/vcfStats.sh ${SAMPLE} ${REF})