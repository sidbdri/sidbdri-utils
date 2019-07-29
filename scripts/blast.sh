#!/usr/bin/env bash

set -o nounset
set -o errexit
#set -o xtrace

function usage {
  cat <<EOT

This script is a wrapper for the ncbi blastn script. It will query the ncbi-blast container with either a sequence or a fasta
file containing multiple reads. The database and index is build by the makembindex command from the blast+ suite.

Usage:
  $(basename $0)
          -q query string or a fasta file contains a number of reads

    Optional:
          -c docker container name which provide the blast service. default: blast
          -f output format. see below for details. see `docker run --rm ncbi/blast blastn -help` for details. default: 0
          -d database name to blast against. default: nt_v5
          -i USEINDEX. default: on

Detail blastn help can be found by running:
> docker run --rm ncbi/blast blastn -help

Useful links:
https://github.com/jarekbryk/localblast
https://github.com/ncbi/docker/tree/master/blast#install-ncbi-provided-blast-databases

Example:
    bash /home/xinhe/Projects/siDBdri-utils/scripts/blast.sh -i -q 'CCCGGAGCCGGTCGCGGCGCACCGCCACGGTGGAAGTGC' -f 0
    bash /home/xinhe/Projects/siDBdri-utils/scripts/blast.sh -i -q 'NTTATATGATAGGGCTCTAAAGATTTTTCAGTATCCAGTACAGTAGCTCTGTAGTAGCTCTGTAGTGCAAATTCATGCACAGTGGAGCTTCAGCCCTTTA' -f '6 qaccver saccver staxid ssciname scomname pident length mismatch gapopen qstart qend sstart send evalue bitscore'
    bash /home/xinhe/Projects/siDBdri-utils/scripts/blast.sh -i -q /home/xinhe/tmp/blast/queries/1202_AS.fasta -f '6 qaccver saccver staxid ssciname scomname pident length mismatch gapopen qstart qend sstart send evalue bitscore'

EOT
}

query=NCAGGTGGGCCTTCCCGGCCGTCCCGGAGCCGGTCGCGGCGCACCGCCACGGTGGAAGTGCGCCCGGCGGCGGCCGGTCGCCGGCCGGGGGGCGGTCCCC
DB=nt_v5
outfmt='0'
USEINDEX=false
QUERY_FOLDER="/srv/data/blast/queries"
TAX_IDS="9606,10090,10116"
DOCKER_CONTAINER=blast

while getopts 'q:c:f:d:i' opt; do
  case ${opt} in
    q)  query=${OPTARG}    ;;
    d)  DB=${OPTARG}        ;;
    f)  outfmt=${OPTARG}     ;;
    c)  DOCKER_CONTAINER=${OPTARG}  ;;
    i)  USEINDEX=true   ;;
    *)  usage; exit         ;;
  esac
done

BASE_DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd -P )



## Check if container exits
[[ ! "$(docker ps -q -f name=${DOCKER_CONTAINER})" ]] && >&2 echo "Error: cannot find running container <${DOCKER_CONTAINER}>" && exit 1


function fastq2fasta {
    local fastq=$1
    local fasta=$2
    local number_of_lines=$3

    head -${number_of_lines} ${fastq} | sed -n '1~4s/^@/>/p;2~4p' > ${fasta}
}

function fastqgz2fasta {
    local fastq=$1
    local fasta=$2
    local number_of_lines=$3
    zcat ${fastq} | head -${number_of_lines} | sed -n '1~4s/^@/>/p;2~4p' > ${fasta}
}


## input is a file, we detect its format for read
if [ -e "${query}" ]; then
#    echo "Input file found: ${query}"
    base_name=$(basename "${query}")
    fmt=${base_name##*.}
    query_name=${base_name%.*}

    tmp="${QUERY_FOLDER}/${query_name}.fasta"

    case ${fmt} in
        fasta)  [[ "${query}" != "${tmp}" ]] && cp ${query} ${tmp};;
#        fastq)  fastq2fasta $query ${tmp} 40 ;;
#        fastq.gz) fastqgz2fasta $query ${tmp} 40 ;;
        *)             echo "Unknown format found: ${fmt}" && exit 1;;
    esac

    query=/blast/queries/${query_name}.fasta

    docker exec ${DOCKER_CONTAINER} /bin/bash -c "blastn -db ${DB} -use_index ${USEINDEX} -outfmt \"${outfmt}\" -query ${query} -taxids ${TAX_IDS}"
else
    docker exec ${DOCKER_CONTAINER} /bin/bash -c "echo ${query} | blastn  -db ${DB} -use_index ${USEINDEX} -outfmt \"${outfmt}\" -taxids ${TAX_IDS}"
fi


