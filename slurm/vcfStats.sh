#!/usr/bin/env bash

#! RUN : sbatch vcfStats.sh <SAMPLE> REF>

#! sbatch directives begin here ###############################
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time 00:30:00
#SBATCH --mail-type=FAIL,END
#SBATCH -p cclake

#SBATCH -o logs/vcfStats_%j.out

. /etc/profile.d/modules.sh                # Leave this line (enables the module command)
module purge                               # Removes all modules still loaded
module load rhel8/default-ccl              # REQUIRED - loads the basic environment

SAMPLE=$1
REF=$2

source ${SAMPLE}.config
STATS_FILE="${SAMPLE}-${REF}.stats"

${BCFTOOLS} stats -s- -f PASS ${SAMPLE}-${REF}.filtered.vcf.gz > output-${REF}.vchk

echo -e "Variants\t"`zcat ${SAMPLE}-${REF}.vcf.gz | grep -v '#' | wc -l` >> ${STATS_FILE}
echo -e "Filtered\t"`zcat ${SAMPLE}-${REF}.filtered.vcf.gz | grep -v '#' | grep -v 'basic' | wc -l` >> ${STATS_FILE}

#echo -e "SNPs\t"`zcat ${SAMPLE}-${REF}.snps.vcf.gz | grep -v '#' | grep -v 'basic' | wc -l` >> ${STATS_FILE}
echo -e "SNPs\t"`grep 'number of SNPs:' output-${REF}.vchk | cut -f 4` >> ${STATS_FILE}
echo -e "snpHom\t"`grep PSC output-${REF}.vchk | tail -1 | cut -f 5` >> ${STATS_FILE}
echo -e "snpHet\t"`grep PSC output-${REF}.vchk | tail -1 | cut -f 6` >> ${STATS_FILE}

#echo -e "InDels\t"`zcat ${SAMPLE}-${REF}.indels.vcf.gz | grep -v '#' | grep -v 'basic' | wc -l` >> ${STATS_FILE}
echo -e "InDels\t"`grep 'number of indels:' output-${REF}.vchk | cut -f 4` >> ${STATS_FILE}
echo -e "insHet\t"`grep PSI output-${REF}.vchk | tail -1 | cut -f 8` >> ${STATS_FILE}
echo -e "delHet\t"`grep PSI output-${REF}.vchk | tail -1 | cut -f 9` >> ${STATS_FILE}
echo -e "insHom\t"`grep PSI output-${REF}.vchk | tail -1 | cut -f 10` >> ${STATS_FILE}
echo -e "delHom\t"`grep PSI output-${REF}.vchk | tail -1 | cut -f 11` >> ${STATS_FILE}


echo -e "dateCompleted\t"`date +"%Y-%m-%d %H:%M:%S"` >> ${SAMPLE}-${REF}.stats

exit;

if [ ! -f ${SAMPLE}-${REF}.filtered.vcf.gz ]; then
	exit 1;
else
  input_size=$(stat -c%s "${SAMPLE}-${REF}.vcf.gz")
  output_size=$(stat -c%s "${SAMPLE}-${REF}.filtered.vcf.gz")
  if [ output_size > input_size ]; then
    #rm -rf ${SAMPLE}-${REF}.g.vcf.gz*
    #rm -rf ${SAMPLE}-${REF}.vcf.gz*
  else
    exit 1;
  fi
fi
