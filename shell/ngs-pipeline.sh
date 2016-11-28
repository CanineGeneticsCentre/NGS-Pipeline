#!/bin/bash

export RESULTS=$HOME/fastq2vcf_runs
export PIPELINE_HOME=$HOME/git/NGS-Pipeline/shell
SAMPLE=$1
STEP=1

if [ "$#" -ge 2 ]; then
	STEP=$2
fi
export START_STEP=${STEP}

if [ ! -f ${RESULTS}/${SAMPLE}.cfg ]; then
	echo "You need to create a configuration file to run the NGS Pipeline for ${SAMPLE} - ${RESULTS}/${SAMPLE}.cfg"
	exit
fi

echo Starting NGS Pipeline for $SAMPLE...
source ${RESULTS}/${SAMPLE}.cfg

if [ ! -e $RESULTS/$SAMPLE ]; then
	mkdir $RESULTS/$SAMPLE
fi

cd $RESULTS/$SAMPLE

if [ ! -e logs ]; then
	mkdir logs
fi

## Convert chunked FASTQ files to aligned SAM files...
##STEP 1
qsub -M ${EMAIL} -V -t 1-${FQ_COUNT} -N ${SAMPLE}.fastq2sam ${PIPELINE_HOME}/step1-fastq2sam.sh ${SAMPLE}

##STEP 2
qsub -M ${EMAIL} -V -t 1-${FQ_COUNT} -hold_jid_ad ${SAMPLE}.fastq2sam -N ${SAMPLE}.sam2bam ${PIPELINE_HOME}/step2-sam2bam.sh ${SAMPLE}

##STEP 3
qsub -M ${EMAIL} -V -hold_jid ${SAMPLE}.sam2bam -N ${SAMPLE}.bam ${PIPELINE_HOME}/step3-bam.sh ${SAMPLE}

##STEP 4
qsub -M ${EMAIL} -V -hold_jid ${SAMPLE}.bam -N ${SAMPLE}.indels ${PIPELINE_HOME}/step4-indels_clean.sh ${SAMPLE}
