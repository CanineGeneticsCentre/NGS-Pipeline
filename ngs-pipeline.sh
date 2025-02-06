#!/usr/bin/env bash

#! RUN : bash ngs-pipeline.sh <SAMPLE> <REF>
#! Eg. : bash ngs-pipeline.sh CS_35365 cf4

SAMPLE=$1
REF=$2
SCRIPTS=`dirname $0`
CFG="${SCRIPTS}/ngs-pipeline-${REF}.config"

[[ -z "$SAMPLE" ]] && { echo "ERROR: No SAMPLE provided for this run"; exit 1; }
[[ -z "$REF" ]] && { echo "ERROR: No REFERENCE provided for this run"; exit 1; }

date +"%d/%m/%Y %H:%M"
echo $SAMPLE $REF;

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
rsync --progress -av --no-links ${WGS}/${SAMPLE}/*.fq.gz ./

# Count how many fastq files we have
COUNT=`ls *.fq.gz | wc -l`
# Work out number of lanes used... assuming pair-end, therefore divide total number of files by 2 - forward and reverse files for each lane
LANES=$((COUNT / 2))

if [ $LANES -le 1 ]; then
  echo ; echo "ERROR - You  need to split the FASTQ files up before running the NGS Pipeline. Please run ***${SCRIPTS}/splitFastq-pipeline.sh ${SAMPLE} ${REF}*** first."; echo
  exit 1;
fi

ls *r_1.fq.gz > fastq.files
if [ ! -s fastq.files ];then
  echo ; echo "ERROR - The FASTQ filenames are not in the required format. Please rename them to end in r_1/r_2.fq.gz"; echo
  exit 1;
fi


mkdir -p $GENOME; 
mv ${SAMPLE}-${REF}.config $GENOME/${SAMPLE}.config
cd $GENOME

scp kcgc:/home/kcgc/tmp/sampleLookups.txt /home/es904/tmp/sampleLookups.txt
grep -i ${SAMPLE}  /home/es904/tmp/sampleLookups.txt | awk '{ print "sample\t" $1 "\nwgsID\t" $3 "\nBioSample\t" $4 }' > ${SAMPLE}-${REF}.stats
echo -e "refGenome\t${REF}" >> ${SAMPLE}-${REF}.stats
echo -e "dateSubmitted\t"`date +"%Y-%m-%d %H:%M:%S"` >> ${SAMPLE}-${REF}.stats

mkdir metrics logs intervals base_recal gvcf

# count raw reads
sbatch -A ${ACCOUNT} -J ${SAMPLE}.countReads ${SCRIPTS}/slurm/countRawReads.sh ${SAMPLE} ${REF}

module load ${GATK}
gatk SplitIntervals -R ${FASTA}/${GENOME}.fasta -L ${INTERVAL_LIST} --scatter-count ${INTERVALS} -O intervals --subdivision-mode BALANCING_WITHOUT_INTERVAL_SUBDIVISION

#BALANCING_WITHOUT_INTERVAL_SUBDIVISION means that --scatter-count is the upper limit for count of intervals
#Need to reassign INTERVALS to the number actually generated
INTERVALS=`ls intervals | wc -l`
sed -i s/INTERVALS=50/INTERVALS=${INTERVALS}/ ${SAMPLE}.config

# fastq2sam
jid1=$(sbatch -A ${ACCOUNT} -J ${SAMPLE}.fastq2sam --array=1-${LANES} ${SCRIPTS}/slurm/make-bam/fastq2sam.sh ${SAMPLE})

# markAdapters
jid2=$(sbatch -A ${ACCOUNT} -J ${SAMPLE}.markAdapters --dependency=afterok:${jid1##* } --array=1-${LANES} ${SCRIPTS}/slurm/make-bam/markAdapters.sh ${SAMPLE})

# sam2fastq - removes ${SAMPLE}.L${LANE}.adaptMarked.bam
jid3=$(sbatch -A ${ACCOUNT} -J ${SAMPLE}.sam2fastq --dependency=afterok:${jid2##* } --array=1-${LANES} ${SCRIPTS}/slurm/make-bam/sam2fastq.sh ${SAMPLE})

# alignfastq
jid4=$(sbatch -A ${ACCOUNT} -J ${SAMPLE}.alignFastq --dependency=afterok:${jid3##* } --array=1-${LANES} ${SCRIPTS}/slurm/make-bam/alignFastq.sh ${SAMPLE})

# mergeBamAlignment (depends on jid1 AND jid4) - removes ${SAMPLE}.L${LANE}.unaligned.bam, ${SAMPLE}.L${LANE}.aligned.bam & ${SAMPLE}.L${LANE}.fastq.gz
jid5=$(sbatch -A ${ACCOUNT} -J ${SAMPLE}.mergeBam --dependency=afterok:${jid4##* } --array=1-${LANES} ${SCRIPTS}/slurm/make-bam/mergeBamAlignment.sh ${SAMPLE})

# markDuplicates - removes lane${LANE}/${SAMPLE}.L${LANE}.merged.bam
# We take advantage of the tool's ability to take multiple BAM inputs and write out a single output
# to avoid having to spend time just merging BAM files.
jid6=$(sbatch -A ${ACCOUNT} -J ${SAMPLE}.markDuplicates --dependency=afterok:${jid5##* } ${SCRIPTS}/slurm/make-bam/markDuplicates.sh ${SAMPLE} ${LANES})

# sortBam - removes ${SAMPLE}.aligned.unsorted.dedup.bam
# Sort BAM file by coordinate order and fix tag values for NM and UQ
jid7=$(sbatch -A ${ACCOUNT} -J ${SAMPLE}.sortBam --dependency=afterok:${jid6##* } ${SCRIPTS}/slurm/make-bam/sortBam.sh ${SAMPLE})

###########################################################
#                   BASE RECALLIBRATION                   #
###########################################################

# Generate Base Quality Score Recalibration (BQSR) model
jid8=$(sbatch -A ${ACCOUNT} -J ${SAMPLE}.BQSR --dependency=afterok:${jid7##* } --array=0-$(($INTERVALS-1)) ${SCRIPTS}/slurm/make-bam/baseRecalibrator.sh ${SAMPLE})

# Merge the recalibration reports resulting from by-interval recalibration
jid9=$(sbatch -A ${ACCOUNT} -J ${SAMPLE}.gatherBQSR --dependency=afterok:${jid8##* } ${SCRIPTS}/slurm/make-bam/gatherBSQR.sh ${SAMPLE})

# Apply the recalibration model by interval
jid10=$(sbatch -A ${ACCOUNT} -J ${SAMPLE}.applyBQSR --dependency=afterok:${jid9##* } --array=0-$(($INTERVALS-1)) ${SCRIPTS}/slurm/make-bam/applyBSQR.sh ${SAMPLE})

# Combine multiple recalibrated BAM files from scattered ApplyRecalibration runs - removes base_recal dir, ${SAMPLE}.sorted.bam, ${SAMPLE}.sorted.bai, ${SAMPLE}.sorted.bam.md5 & ${SAMPLE}.bsqr.out
jid11=$(sbatch -A ${ACCOUNT} -J ${SAMPLE}.BAM --dependency=afterok:${jid10##* } ${SCRIPTS}/slurm/make-bam/bqsrMerge.sh ${SAMPLE} ${REF})


###########################################################
#                         METRICS                         #
###########################################################

# QC - Yield metrics
tmp=$(sbatch -A ${ACCOUNT} -J ${SAMPLE}.yieldQC --dependency=afterok:${jid1##* } --array=1-${LANES} ${SCRIPTS}/slurm/make-bam/qcYieldMetrics.sh ${SAMPLE})

# QC - Lane/RG metrics
tmp=$(sbatch -A ${ACCOUNT} -J ${SAMPLE}.laneQC --dependency=afterok:${jid5##* } --array=1-${LANES} ${SCRIPTS}/slurm/make-bam/qcLaneMetrics.sh ${SAMPLE})

# QC - Final BAM file metrics - insert size, flagstats etc.
jid12=$(sbatch -A ${ACCOUNT} -J ${SAMPLE}.bamQC --dependency=afterok:${jid11##* } ${SCRIPTS}/slurm/make-bam/qcBamMetrics.sh ${SAMPLE} ${REF})


###########################################################
#                        gVCF File                        #
###########################################################

# Create gvcf files with HaplotypeCaller
jid13=$(sbatch -A ${ACCOUNT} -J ${SAMPLE}.HC --dependency=afterok:${jid11##* } --array=0-$(($INTERVALS-1)) ${SCRIPTS}/slurm/make-gvcf/haplotypeCaller.sh ${SAMPLE} ${REF} ${PCR_MODEL})

# Merge gVCF files into single gVCF
#jid14=$(sbatch -A ${ACCOUNT} -J ${SAMPLE}.GVCF --dependency=afterok:${jid13##* } ${SCRIPTS}/slurm/make-gvcf/combineGvcf.sh ${SAMPLE} ${INTERVALS} ${REF})
jid14=$(sbatch -A ${ACCOUNT} -J ${SAMPLE}.GVCF --dependency=afterok:${jid13##* } ${SCRIPTS}/slurm/make-gvcf/mergeGvcfs.sh ${SAMPLE})


###########################################################
#                     VCF Files/Stats                     #
###########################################################

# Convert gVCF file to VCF
jid15=$(sbatch -A ${ACCOUNT} -J ${SAMPLE}.VCF --dependency=afterok:${jid14##* } ${SCRIPTS}/slurm/gvcf2vcf.sh ${SAMPLE} ${REF})

# Filter VCF and split into snps/indels
jid16=$(sbatch -A ${ACCOUNT} -J ${SAMPLE}.filterVCF --dependency=afterok:${jid15##* } ${SCRIPTS}/slurm/filterVcf.sh ${SAMPLE} ${REF})

# Output VCF stats
tmp=$(sbatch -A ${ACCOUNT} -J ${SAMPLE}.VCFstats --dependency=afterok:${jid16##* } ${SCRIPTS}/slurm/vcfStats.sh ${SAMPLE} ${REF})



###########################################################
#                         OUTPUT                          #
###########################################################
echo "BAM file = job-${jid11##* }"
echo "gVCF file = job-${jid14##* }"
echo "VCF file = job-${jid16##* }"
echo "bamMetrics-${jid12##* }"
