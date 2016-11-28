#!/bin/bash -l

#$ -b y
#$ -cwd
#$ -l h_vmem=3G
#$ -pe smp 6
#$ -o logs/step4.out
#$ -e logs/step4.err
#$ -m ae

if [ $START_STEP -gt 4  ]; then
	exit
fi

module load apps/picard

SAMPLE=$1
source ${RESULTS}/${SAMPLE}.cfg

##RealignerTargetCreator
gatk -T RealignerTargetCreator -nt 6 -R ${CF3}/canfam3.fasta -I ${SAMPLE}.bam -o ${SAMPLE}.bam.intervals -mismatch 0.0 -S LENIENT --filter_mismatching_base_and_quals

##IndelRealigner
gatk -T IndelRealigner -R ${CF3}/canfam3.fasta -I ${SAMPLE}.bam -targetIntervals ${SAMPLE}.bam.intervals -o ${SAMPLE}.clean.bam -S LENIENT --filter_mismatching_base_and_quals

### if size of the clean bam file >= 10% of bam then delete file
bam_size=$(wc -c < ${RESULTS}/${SAMPLE}/${SAMPLE}.bam)
bam_size1=$(wc -c < ${RESULTS}/${SAMPLE}/${SAMPLE}.clean.bam)
if [ $bam_size1 -ge $(( ${bam_size}/10 )) ]; then
	echo "Deleting ${RESULTS}/${SAMPLE}/${SAMPLE}.ba*"
	rm -rf ${RESULTS}/${SAMPLE}/${SAMPLE}.bam
	rm -rf ${RESULTS}/${SAMPLE}/${SAMPLE}.bai
fi


##BaseRecalibrator
gatk -T BaseRecalibrator -nct 6 -R ${CF3}/canfam3.fasta -I ${SAMPLE}.clean.bam -knownSites ${CF3}/canfam3_all_snps.vcf -o ${SAMPLE}.recal.grp -S LENIENT

##PrintReads
gatk -T PrintReads -nct 6 -R ${CF3}/canfam3.fasta -I ${SAMPLE}.clean.bam -BQSR ${SAMPLE}.recal.grp -o ${SAMPLE}.final.bam -S LENIENT

### if size of the final bam file >= 10% of clean bam then delete file
bam_size2=$(wc -c < ${RESULTS}/${SAMPLE}/${SAMPLE}.final.bam)
if [ $bam_size1 -ge $(( ${bam_size1}/10 )) ]; then
	echo "Deleting ${RESULTS}/${SAMPLE}/${SAMPLE}.clean.ba*"
	rm -rf ${RESULTS}/${SAMPLE}/${SAMPLE}.clean.bam
	rm -rf ${RESULTS}/${SAMPLE}/${SAMPLE}.clean.bai
fi
