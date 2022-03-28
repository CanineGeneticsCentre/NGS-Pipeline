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

mkdir -p $SAMPLE/logs; cd $SAMPLE
cp $CFG $SAMPLE.config; source $SAMPLE.config

#PCR_FREE=false;
PCR_MODEL='CONSERVATIVE';

printf "Is the data PCR free?\n"
printf "\t1. No [default]\n"
printf "\t2. Yes\n"
# Assign input value into a variable
read answer

if [[ -v $answer && $answer == "2" ]]; then
    #PCR_FREE=true;
    PCR_MODEL='NONE';
fi

# Copy files from RCS... will only copy if files no present or rcs version is newer
rsync --progress -av ${WGS}/${SAMPLE}/*.fq.gz ./
#rsync --progress -av ${WGS}/${SAMPLE}/*.fastq.gz ./

# Count how many fastq files we have
COUNT=`ls *.fq.gz | wc -l`
# Work out number of lanes used... assuming pair-end, therefore divide total number of files by 2 - forward and reverse files for each lane
LANES=$((COUNT / 2))

if [ $LANES -le 1 ]; then
  echo ; echo "ERROR - You  need to split the FASTQ files up before running the NGS Pipeline. Please run ***${SCRIPTS}/splitFastq-pipeline.sh ${SAMPLE} ${REF}*** first."; echo
  exit 1;
fi


mkdir -p $GENOME; 
mv $SAMPLE.config $GENOME/
cd $GENOME

# generate seqeunce groups for future scatter/gather steps.
perl ${SCRIPTS}/perl/createSeqGroups.pl ${DICT}
INTERVALS=`wc -l sequence_grouping.txt | awk '{print $1}'`

BARCODE=`ls *.fq.gz | head -1 | xargs -n 1 zcat 2>/dev/null | head -1 | cut -d':' -f 1-3`

# fastq2bam - Submit job array to align samples to ref genome
jid1=$(sbatch -A ${ACCOUNT} -J ${SAMPLE}.fastq2bam --array=1-${LANES} ${SCRIPTS}/slurm/fastq2bam.sh ${SAMPLE})

# addRGinfo - add Read Group information to aligned BAM files
jid2=$(sbatch -A ${ACCOUNT} -J ${SAMPLE}.rg --dependency=afterok:${jid1##* } --array=1-${LANES} ${SCRIPTS}/slurm/addRGinfo.sh ${SAMPLE} ${BARCODE})

# mark duplicates
# We take advantage of the tool's ability to take multiple BAM inputs and write out a single output
# to avoid having to spend time just merging BAM files.
jid3=$(sbatch -A ${ACCOUNT} -J ${SAMPLE}.markDuplicates --dependency=afterok:${jid2##* } ${SCRIPTS}/slurm/markDuplicates.sh ${SAMPLE} ${LANES})

# Sort BAM file by coordinate order and fix tag values for NM and UQ
jid4=$(sbatch -A ${ACCOUNT} -J ${SAMPLE}.sortSam --dependency=afterok:${jid3##* } ${SCRIPTS}/slurm/sortSam.sh ${SAMPLE})

# Generate Base Quality Score Recalibration (BQSR) model
jid5=$(sbatch -A ${ACCOUNT} -J ${SAMPLE}.BQSR --dependency=afterok:${jid4##* } --array=1-${INTERVALS} ${SCRIPTS}/slurm/baseRecalibrator.sh ${SAMPLE})

# Merge the recalibration reports resulting from by-interval recalibration
jid6=$(sbatch -A ${ACCOUNT} -J ${SAMPLE}.gatherBQSR --dependency=afterok:${jid5##* } ${SCRIPTS}/slurm/gatherBSQR.sh ${SAMPLE})


INTERVALS=`wc -l sequence_grouping_with_unmapped.txt | awk '{print $1}'`
# Apply the recalibration model by interval
jid7=$(sbatch -A ${ACCOUNT} -J ${SAMPLE}.applyBQSR --dependency=afterok:${jid6##* } --array=1-${INTERVALS} ${SCRIPTS}/slurm/applyBSQR.sh ${SAMPLE})

# Combine multiple recalibrated BAM files from scattered ApplyRecalibration runs
jid8=$(sbatch -A ${ACCOUNT} -J ${SAMPLE}.BAM --dependency=afterok:${jid7##* } ${SCRIPTS}/slurm/mergeBam.sh ${SAMPLE} ${INTERVALS} ${REF})


INTERVALS=`wc -l sequence_grouping.txt | awk '{print $1}'`
# Create gvcf files with HaplotypeCaller
jid9=$(sbatch -A ${ACCOUNT} -J ${SAMPLE}.HC --dependency=afterok:${jid8##* } --array=1-${INTERVALS} ${SCRIPTS}/slurm/haplotypeCaller.sh ${SAMPLE} ${REF} ${PCR_MODEL})

# Merge gVCF files into single gVCF
jid10=$(sbatch -A ${ACCOUNT} -J ${SAMPLE}.GVCF --dependency=afterok:${jid9##* } ${SCRIPTS}/slurm/combineGvcf.sh ${SAMPLE} ${INTERVALS} ${REF})



# Generate metrics for final BAM file - insert size, flagstats etc.
jid11=$(sbatch -A ${ACCOUNT} -J ${SAMPLE}.metrics --dependency=afterok:${jid8##* } ${SCRIPTS}/slurm/bamMetrics.sh ${SAMPLE} ${REF})



#jid3=$(sbatch -J ${SAMPLE}.validateSam --dependency=afterok:${jid2##* } ${SCRIPTS}/slurm/validateSam.sh ${SAMPLE})
#
#jid5=$(sbatch -J ${SAMPLE}.realignReads --dependency=afterok:${jid4##* } ${SCRIPTS}/slurm/realignReads.sh ${SAMPLE})
#jid6=$(sbatch -J ${SAMPLE}.realignIndels --dependency=afterok:${jid5##* } ${SCRIPTS}/slurm/realignIndels.sh ${SAMPLE})
#jid7=$(sbatch -J ${SAMPLE}.baseRecalibrator --dependency=afterok:${jid6##* } ${SCRIPTS}/slurm/baseRecalibrator.sh ${SAMPLE})
#jid8=$(sbatch -J ${SAMPLE}.printReads --dependency=afterok:${jid7##* } ${SCRIPTS}/slurm/printReads.sh ${SAMPLE})
#jid9=$(sbatch -J ${SAMPLE}.CollectMetrics --dependency=afterok:${jid8##* } ${SCRIPTS}/slurm/collectMetrics.sh ${SAMPLE})
#jid10=$(sbatch -J ${SAMPLE}.validateSam2 --dependency=afterok:${jid8##* } ${SCRIPTS}/slurm/validateSam2.sh ${SAMPLE})
#jid11=$(sbatch -J ${SAMPLE}.FinalStep --dependency=afterok:${jid9##* }:${jid10##* } ${SCRIPTS}/slurm/finalStep.sh ${SAMPLE})

echo "BAM file = job-${jid8##* }"
echo "gVCF file = job-${jid10##* }"
echo "Metrics = job-${jid11##* }"
