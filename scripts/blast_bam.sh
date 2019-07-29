#!/usr/bin/env bash

set -o nounset
set -o errexit
set -o xtrace

function cleanup {
#   echo "Killing all sub-processes..."
   kill -- -$$
}

trap exit INT
trap cleanup EXIT

export BASE_DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd -P )
#source "${BASE_DIR}/../includes.sh"

function usage {
  cat <<EOT

Usage:
  $(basename $0)
        -b bam file to be process
        -r randomly select a number of reads from a range in the bam file for blast
            -g the range from which to select the random reads
            -n the number of random reads to be selected
           * Incompatible with:  -i

        -i included reads. example: '1-3,4,5-6' will select reads  1 to 6 from bam file for blast
           * Incompatible with:  -r

   Optional parameter:
        -q queries directory. This is the directory mounted to the blast docker container containing quesies of
            fasta files. default /srv/data/blast/queries
        -t temp directory for intermediate files. default: /tmp/blast_tmp
        -o blast result file name. default: -t/-b.txt

Example:
    # select read 1 and read 2 from bam for blast
    bash blast_bam.sh -b /srv/data/blast/fasta/1202_As___rat___filtered.bam -i 1,2
    bash blast_bam.sh -b /srv/data/blast/fasta/1202_As___rat___filtered.bam -i 1,2 -o  /tmp/blast_tmp/out.txt

    # select read 3 to read 5 from bam for blast
    bash blast_bam.sh -b /srv/data/blast/fasta/1202_As___rat___filtered.bam -i 3-5

    # select 2 random reads from read 1 to read 10 for blast
    bash blast_bam.sh -b /srv/data/blast/fasta/1202_As___rat___filtered.bam -r -g 1-10 -n 2

EOT
}

## default value
QUERY_DIR=/srv/data/blast/queries

BAM_FILE=/srv/data/blast/fasta/1202_As___rat___filtered.bam
TMP_DIR=/tmp/blast_tmp/


OUT=${TMP_DIR}/`basename ${BAM_FILE} | sed 's/.bam/.blast/' `

fasta_file=${TMP_DIR}/`basename ${BAM_FILE} | sed 's/.bam/.fasta/' `
sed_pattern_file=${TMP_DIR}/`basename ${BAM_FILE} | sed 's/.bam/.sed/' `
query_file=${QUERY_DIR}/`basename ${BAM_FILE} | sed 's/.bam/.fasta/'`

RDM="false"
RDM_RANGE='1-5'
RDM_NUM_READ=3
INCLUDE_RECORD='1-5,6,7,9-10'


## read parameter
while getopts 'b:q:t:o:rg:n:i:' opt; do
  case ${opt} in
    b)  BAM_FILE=${OPTARG}    ;;
    q)  QUERY_DIR=${OPTARG}    ;;
    t)  TMP_DIR=${OPTARG}    ;;
    o)  OUT=${OPTARG}        ;;
    r)  RDM="true"            ;;
    g)  RDM_RANGE=${OPTARG}            ;;
    n)  RDM_NUM_READ=${OPTARG}            ;;
    i)  INCLUDE_RECORD=${OPTARG}            ;;
    *)  usage; exit         ;;
  esac
done



mkdir -p ${TMP_DIR}

[[ ! -e ${BAM_FILE} ]] && echo 'Error: ${BAM_FILE} not found.' && exit
##we generate a tmp fasta file contains all reads from bam

[[ ! -e ${fasta_file} ]] && samtools bam2fq ${BAM_FILE} | sed -n '1~4s/^@/>/p;2~4p' > ${fasta_file}



if [ "${RDM}" = "true" ]; then
    ## we randomly select a number of read from a range
    echo "Selecting reads from fasta to create query file..."

    records=`shuf -i ${RDM_RANGE} -n ${RDM_NUM_READ} | sort -n | uniq`
else
    ## we select reads specfied by INCLUDE
    non_range=`echo ${INCLUDE_RECORD} | awk 'BEGIN{RS=","} {print}' | grep -v '-'`
    range=`echo "${non_range}"|head -1`
    if [ -n "$(echo ${INCLUDE_RECORD} | grep '-')" ]; then
        ##there is range
        range=`echo ${INCLUDE_RECORD} | awk 'BEGIN{RS=","} {print}' | grep '-' | sed 's/-/ 1 /' | xargs -L 1 seq`
    fi
    records=`printf "$range\n$non_range" | sort -n | uniq`
fi

max_row=`echo "${records}" | tail -1 | awk '{print 2*$0+1}'`
pattern=`echo "${records}" | awk '{print 2*$0-1"p;"2*$0"p;" }'|  paste -sd' '`

echo ${pattern}${max_row}q > ${sed_pattern_file}
sed -n -f ${sed_pattern_file} ${fasta_file} > ${query_file}

cat ${query_file}

#/usr/bin/time -f "%e"
bash ${BASE_DIR}/blast.sh -i -q ${query_file} -f '6 qaccver saccver staxid ssciname pident evalue bitscore' > ${OUT}
echo "blast results saved in ${OUT}"

## Generate a summary
Rscript -e "
library(magrittr,quietly = T,warn.conflicts = F); \
library(dplyr,quietly = T,warn.conflicts = F); \
read.table(file = '${OUT}',sep = '\t') %>% \
    group_by(V1) %>% summarise(s=head(V4,1)) %>% \
    pull(s) %>% table() %>% prop.table()"