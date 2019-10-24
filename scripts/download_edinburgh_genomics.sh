#!/usr/bin/env bash


set -o nounset
set -o errexit
#set -o xtrace

function cleanup {
#   echo "Killing all sub-processes..."
   kill -- -$$
}

trap exit INT
trap cleanup EXIT

export BASE_DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd -P )

function usage {
  cat <<EOT

Usage:
  $(basename $0)
        -u Edinburgh Genomics username
        -p Edinburgh Genomics password
        -d project path
        -o download folder
Example:
    bash $(basename $0) -u xin -p 123455 -d "11907_Hardingham_Giles/raw_data/20191011" -o /srv/data/ghardingham/contralateral_cortex_zoeb

EOT
}

## default value
#USERNAME=""
#PASSWORD=""
#PROJECT_PATH="11907_Hardingham_Giles/raw_data/20191011"
#OUT_DIR="/srv/data/ghardingham/contralateral_cortex_zoeb"

## read parameter
while getopts 'u:p:d:o:' opt; do
  case ${opt} in
    u)  USERNAME=${OPTARG}    ;;
    p)  PASSWORD=${OPTARG}    ;;
    d)  PROJECT_PATH=${OPTARG}    ;;
    o)  OUT_DIR=${OPTARG}    ;;
    *)  usage; exit         ;;
  esac
done


mkdir -p ${OUT_DIR}
cd ${OUT_DIR}
export ASPERA_SCP_PASS=${PASSWORD}
nohup ascp -P 33001 -O 33001 ${USERNAME}@transfer.genomics.ed.ac.uk:${PROJECT_PATH} . > ${OUT_DIR}/download_log.txt &
wait

echo "done!" >> ${OUT_DIR}/download_log.txt