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


    This script turns all the .genozip file in a specific folder into .fastq.gz format.
    Here is a brief summary of how it works:
    For a given folder under contains all the sample folders, the script
        1. check if every .genozip files in the given folder of a corresponding .fastq.gz file
        2. uncompress those .genozip files that don't have corresponding .fastq.gz files


    Given how the data is structured on sidb, here is a example workflow:
        >cd /srv/data/fastq_temp
        # we copy the project folder containing all the sample .fastq.gz files
        >mkdir -p nrf2_ich_jloan
        >cd /srv/data/fastq_temp/nrf2_ich_jloan
        >find /srv/data/ghardingham/nrf2_ich_jloan -type d | tail -n +2  | sort | head -2 | xargs -P 31 -n 1 -t -I % bash -c "cp -R % ."
        >bash /home/xinhe/Projects/sidbdri-utils/scripts/prepare_fastq_from_genozip.sh -i /srv/data/fastq_temp/nrf2_ich_jloan -f

    Here is an example output of the script with only two samples of the nrf2_ich_jloan project:
        >bash /home/xinhe/Projects/sidbdri-utils/scripts/prepare_fastq_from_genozip.sh -i /srv/data/fastq_temp/nrf2_ich_jloan -f
            Checking if every .genozip file has the corresponding .fastq.gz file...
            /srv/data/fastq_temp/nrf2_ich_jloan/10_Flox_Ma_ICH/10_Flox_Ma_ICH_R1_merged.fastq.genozip... failed!
            /srv/data/fastq_temp/nrf2_ich_jloan/10_Flox_Mg_ICH/10_Flox_Mg_ICH_R1_merged.fastq.genozip... pass!
            Script will uncompress the following .genozip files in /srv/data/fastq_temp/nrf2_ich_jloan...
            /srv/data/fastq_temp/nrf2_ich_jloan/10_Flox_Ma_ICH/10_Flox_Ma_ICH_R1_merged.fastq.genozip
            uncompressing with genounzip...
            bash -c uncompress_genozip /srv/data/fastq_temp/nrf2_ich_jloan/10_Flox_Ma_ICH/10_Flox_Ma_ICH_R1_merged.fastq.genozip true
            Checking if every .genozip file has the corresponding .fastq.gz file...
            /srv/data/fastq_temp/nrf2_ich_jloan/10_Flox_Ma_ICH/10_Flox_Ma_ICH_R1_merged.fastq.genozip... pass!
            /srv/data/fastq_temp/nrf2_ich_jloan/10_Flox_Mg_ICH/10_Flox_Mg_ICH_R1_merged.fastq.genozip... pass!
            Script finished.


    Usage:
      $(basename $0)
            -d folder which contains the sample folders
            -f genozip will overwrite existing genozip file without warning
            -p number of cores. default: 16
            -h see help page
    Example:
        bash /home/xinhe/Projects/sidbdri-utils/scripts/prepare_fastq_from_genozip.sh -i /srv/data/fastq_temp/nrf2_ich_jloan -f -p 10

EOT
}

set -a
function uncompress_genozip {
    local file=$1
    local force=${2:-false}

    local para=''
    [[ "${force}" == "true" ]] && para='-f'

    para=${para}" --quiet"

    if [  -f "${file}" ]; then
        # after compress with genozip, we uncompress it to check its md5 hash to make sure it is the same as the original one
        genounzip ${para} ${file}
    else
        echo "genozip file <${file}> does not exist."
        return 1
    fi
}
set +a

## default value
SAMPLE_DIR='/srv/data/fastq_temp'
NUM_CORES=16
OVERWRITE_EXISTING=false

## read parameter
while getopts 'i:p:f' opt; do
  case ${opt} in
    i)  SAMPLE_DIR=${OPTARG}    ;;
    p)  NUM_CORES=${OPTARG}     ;;
    f)  OVERWRITE_EXISTING="true"   ;;
    *)  usage; exit         ;;
  esac
done


function confirm {
    local msg=${1:-"Please double check. Do you wish to continue?"}
    echo ${msg}
    select yn in "Yes" "No"; do
        case $yn in
            Yes ) echo 'yes'; break;;
            No ) echo 'script terminated!'; exit;;
            *) echo "script terminated!"; exit;;
        esac
    done
}



[[ ! "${SAMPLE_DIR}" =~ ^/srv/data/.* ]] && echo "Error: GENOZIP_DATA_DIR can only be under /srv/data/" && exit 1
[[ ! -d "${SAMPLE_DIR}"  ]] && echo "Error: SAMPLE_DIR <${SAMPLE_DIR}> does not exist. " && usage && exit 1


FILES_TO_UNCOMPRESS=()
echo "Checking if every .genozip file has the corresponding .fastq.gz file..."
#echo "The following files will be checked:"
#echo "--------------------------"
#find ${SAMPLE_DIR} -type f -name "*.genozip" | sort
#echo "--------------------------"
for file in `find ${SAMPLE_DIR} -type f -name "*.genozip" | sort`; do
    if [ -f "${file%.genozip}.gz" ]; then
        echo "${file}... pass!"
#        FILES_TO_UNCOMPRESS+=(${file})
    else
        echo "${file}... failed!"
        FILES_TO_UNCOMPRESS+=(${file})
    fi
done

if [ "${#FILES_TO_UNCOMPRESS[@]}" == '0' ]; then
    echo "Every .genozip file has the corresponding .fastq.gz file, nothing needs to be done."
    exit 0
fi


echo "Script will uncompress the following .genozip files in ${SAMPLE_DIR}..."
for file in ${!FILES_TO_UNCOMPRESS[@]}; do
    echo ${FILES_TO_UNCOMPRESS[$file]}
done
#confirm

echo "uncompressing with genounzip..."
for file in ${!FILES_TO_UNCOMPRESS[@]}; do
    echo "uncompress_genozip ${FILES_TO_UNCOMPRESS[$file]} ${OVERWRITE_EXISTING}"
done | xargs -t -n 1 -P ${NUM_CORES} -I % bash -c "%"


echo "Checking if every .genozip file has the corresponding .fastq.gz file..."
for file in `find ${SAMPLE_DIR} -type f -name "*.genozip" | sort`; do
    [[ -f "${file%.genozip}.gz" ]] && echo "${file}... pass!" || echo "${file}... failed!"
done
echo "Script finished."

#find /srv/data/ghardingham/nrf2_ich_jloan -type d | tail -n +2 | xargs -P 31 -n 1 -t -I % bash -c "cp -R % ."