#!/usr/bin/env bash

#! RUN : bash ngs-pipeline.sh <SAMPLE> <REF>
#! Eg. : bash ngs-pipeline.sh CS_35365 cf4

SAMPLE=$1
REF=$2
SCRIPTS=`dirname $0`
CFG="${SCRIPTS}/ngs-pipeline-${REF}.config"

[[ -z "$SAMPLE" ]] && { echo "ERROR: No SAMPLE provided for this run"; exit 1; }
[[ -z "$REF" ]] && { echo "ERROR: No REFERENCE provided for this run"; exit 1; }

date

PCR_MODEL='CONSERVATIVE';

printf "Is the data PCR free?\n"
printf "\t1. No [default]\n"
printf "\t2. Yes\n"
# Assign input value into a variable
read answer

if [[ -n $answer && $answer == "2" ]]; then
    PCR_MODEL='NONE';
fi

mkdir -p $SAMPLE; cd $SAMPLE
cp $CFG ${SAMPLE}-${REF}.config; source ${SAMPLE}-${REF}.config

# Copy files from RCS... will only copy if files no present or rcs version is newer
rsync --progress -av ${WGS}/${SAMPLE}/*.fq.gz ./

# Count how many fastq files we have
COUNT=`ls *.fq.gz | wc -l`
# Work out number of lanes used... assuming pair-end, therefore divide total number of files by 2 - forward and reverse files for each lane
LANES=$((COUNT / 2))

if [ $LANES -le 1 ]; then
  echo ; echo "ERROR - You  need to split the FASTQ files up before running the NGS Pipeline. Please run ***${SCRIPTS}/splitFastq-pipeline.sh ${SAMPLE} ${REF}*** first."; echo
  exit 1;
fi

#BARCODE=`ls *.fq.gz | head -1 | xargs -n 1 zcat 2>/dev/null | head -1 | cut -d':' -f 1-3`

mkdir -p $GENOME; 
mv ${SAMPLE}-${REF}.config $GENOME/${SAMPLE}.config
cd $GENOME

mkdir metrics logs tmp_files intervals

module load ${GATK}
gatk SplitIntervals -R ${FASTA}/${GENOME}.fasta -L ${INTERVAL_LIST} --scatter-count ${INTERVALS} -O intervals --subdivision-mode BALANCING_WITHOUT_INTERVAL_SUBDIVISION

# fastq2sam
jid1=$(sbatch -A ${ACCOUNT} -J ${SAMPLE}.fastq2sam --array=1-${LANES} ${SCRIPTS}/slurm/make-bam/fastq2sam.sh ${SAMPLE})

# QC - Yield metrics
tmp=$(sbatch -A ${ACCOUNT} -J ${SAMPLE}.yield --dependency=afterok:${jid1##* } --array=1-${LANES} ${SCRIPTS}/slurm/make-bam/qcYieldMetrics.sh ${SAMPLE})

# markAdapters
jid2=$(sbatch -A ${ACCOUNT} -J ${SAMPLE}.markAdapters --dependency=afterok:${jid1##* } --array=1-${LANES} ${SCRIPTS}/slurm/make-bam/markAdapters.sh ${SAMPLE})

# sam2fastq
jid3=$(sbatch -A ${ACCOUNT} -J ${SAMPLE}.sam2fastq --dependency=afterok:${jid2##* } --array=1-${LANES} ${SCRIPTS}/slurm/make-bam/sam2fastq.sh ${SAMPLE})

# alignfastq (depends on jid1 AND jid3)
jid4=$(sbatch -A ${ACCOUNT} -J ${SAMPLE}.alignFastq --dependency=afterok:${jid3##* } --array=1-${LANES} ${SCRIPTS}/slurm/make-bam/alignFastq.sh ${SAMPLE})
