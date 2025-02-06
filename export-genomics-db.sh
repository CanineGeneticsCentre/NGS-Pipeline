#!/usr/bin/env bash

#! RUN : bash export-genomics-db.sh <REF> [<CHR>]
#! Eg. : bash export-genomics-db.sh cf4 [38]

REF=$1
CHR=$2
SCRIPTS=`dirname $0`
CFG="${SCRIPTS}/ngs-pipeline-${REF}.config"

[[ -z "$REF" ]] && { echo "ERROR: No REFERENCE provided for this run"; exit 1; }

source ${CFG};
mkdir -p ${GENOME}/logs; cd $GENOME
mkdir -p snpEff
cp $CFG ${REF}.config; source ${REF}.config

# Need to submit a separate job array per chromosome to make jobs small enough to run and complete
# If chromosome is provided then just submit that job, else loop round all chrs in the genomicsDB.intervals file
# Currently will not work for Y/MT/Un on their own...
if [ "$CHR" ]; then
  if [[ ${#CHR} -lt 4 ]] ; then
    CHR="chr"${CHR}
  fi
  LINE=$(grep -n "^${CHR}$" ${FASTA}/${REF}-genomicsDB.intervals | cut -f 1 -d":")
  SEQ="${LINE} ${LINE}"
else
  INTERVALS=`wc -l ${FASTA}/${REF}-genomicsDB.intervals | awk '{print $1}'`
  SEQ="1 ${INTERVALS}"
fi

for i in `seq ${SEQ}`; do
  CHR=$(head -$i ${FASTA}/${REF}-genomicsDB.intervals | tail -1 | cut -f 1 -d' ' | cut -d'_' -f 1 | cut -f 1 -d':')
  if [[ ${#CHR} -lt 4 ]] ; then
    CHR="chr"${CHR}
  fi

  # loop over each occurence of CHR in the dictionary and chunk as necessary
  rm -rf ${CHR}.intervals; touch ${CHR}.intervals
  grep $CHR ${DICT} | while read -r sq SN LN line ; do 
    LEN=$(echo $LN | cut -f 2 -d":");
    SN=$(echo $SN | cut -f 2 -d":");
    if [ `echo $SN | cut -f 1 -d '_'` == ${CHR} ]; then
      perl ${SCRIPTS}/perl/chunk-chr.pl ${SN} ${LEN} >> ${CHR}.intervals
    fi
  done

  INTERVALS=`wc -l ${CHR}.intervals | awk '{print $1}'`
  ARRAY="1-${INTERVALS}"

  mkdir ${CHR}; touch ${CHR}/files.list
  for j in `seq 1 ${INTERVALS}`; do echo ${CHR}/${REF}-${CHR}-$j.filtered.vcf.gz >> ${CHR}/files.list; done

  #EXPORT from genomicsDB per chr/chunk
  jid=$(sbatch -A ${ACCOUNT} -J ${REF}-${CHR}.VCF --array=${ARRAY} --export=SCRIPTS=${SCRIPTS} ${SCRIPTS}/slurm/make-vcf/GenomicsDB2vcf.sh ${REF} ${CHR})
  
  #FILTER each chr-chunk file
  jid2=$(sbatch -A ${ACCOUNT} -J ${REF}-${CHR}.filterVCF --array=${ARRAY} --dependency=afterok:${jid##* } ${SCRIPTS}/slurm/make-vcf/filterMultiVcf.sh ${REF} ${CHR})

  #CONCAT each chr-chunk file into a single chr VCF file
  jid3=$(sbatch -A ${ACCOUNT} -J ${REF}-${CHR}.concatVCF --dependency=afterok:${jid2##* } ${SCRIPTS}/slurm/make-vcf/concatVcf.sh ${REF} ${CHR})

  #ANNOTATE echo chr
  sbatch -A ${ACCOUNT} -J ${REF}-${CHR}.snpEff --dependency=afterok:${jid3##* } ${SCRIPTS}/slurm/make-vcf/annotateVcf.sh ${REF} ${CHR};
done
exit;





# If chromosome is NOT set then use all lines in intervals list to scatter the update job
if [ -z "$CHR" ]; then
  INTERVALS=`wc -l ${FASTA}/${REF}-genomicsDB.intervals | awk '{print $1}'`
  ARRAY="1-${INTERVALS}"
else
  MIN=`grep -n "^$CHR" $FASTA/${REF}-genomicsDB.intervals | awk  -F':' ' { print $1 } ' | head -1`
  MAX=`grep -n "^$CHR" $FASTA/${REF}-genomicsDB.intervals | awk  -F':' ' { print $1 } ' | tail -1`
  ARRAY="${MIN}-${MAX}"
fi

#echo sbatch -A ${ACCOUNT} -J ${REF}.VCF --array=${ARRAY} --export=SCRIPTS=${SCRIPTS} ${SCRIPTS}/slurm/GenomicsDB2vcf.sh ${REF}
jid=$(sbatch -A ${ACCOUNT} -J ${REF}.VCF --array=${ARRAY} --export=SCRIPTS=${SCRIPTS},REF=${REF} ${SCRIPTS}/slurm/GenomicsDB2vcf.sh)

jid2=$(sbatch -A ${ACCOUNT} -J ${REF}.filterVCF --array=${ARRAY} --dependency=afterok:${jid##* } ${SCRIPTS}/slurm/filterMultiVcf.sh ${REF})

if [ -z "$CHR" ]; then
  for chr in `cut -f 1 -d':' ${FASTA}/${REF}-genomicsDB.intervals | cut -f 1 -d' ' | cut -d'_' -f 1 | sort -n | uniq`; do
    sbatch -A ${ACCOUNT} -J ${REF}-$chr.snpEff --dependency=afterok:${jid2##* } --export=SCRIPTS=${SCRIPTS},REF=${REF} ${SCRIPTS}/slurm/annotateVcf.sh $chr;
  done;
else
  sbatch -A ${ACCOUNT} -J ${REF}-${CHR}.snpEff --dependency=afterok:${jid2##* } --export=SCRIPTS=${SCRIPTS},REF=${REF} ${SCRIPTS}/slurm/annotateVcf.sh $CHR;
fi

