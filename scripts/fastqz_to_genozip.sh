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


    This script turns all the .fastq.gz file in a specific folder into .genozip format.
    Here is a brief summary of how it works:
    For a given folder containing all the sample folders, the script
        1. uncompress all the .fastq.gz files into .fastq format
        2. compress the .fastq.gz into .genozip
        3. uncompress .genozip into .fastq format
        4. md5sum all the .fastq files from step 1 and 3
        5. check if every md5sum appear exactly twice
        6. delete all the fastq.gz files
        7. delete all the fastq files
        8. change the owner group of all files in the give folder to sidbdri

    This script is to be run manually as it requires user input to confirm the action in step 1,6 and 7.
    User should take EXTRA CARE when confirming at step 6 and 7 to make sure the script is deleting the correct files.
    User should MAKE SURE the script run successfully to the end. The compress/uncompress process is not 'atomic'.
    Thus the script terminates when a file is being compress/uncompress, it will leave a half complete file in the folder.


    Given how the data is structured on sidb, here is a example workflow:
        >cd /srv/data/fastq_temp
        # we copy the project folder containing all the sample .fastq.gz files as test
        >cp -R /srv/data/ghardingham/nrf2_ich_jloan /srv/data/fastq_temp/
        >cd /srv/data/fastq_temp/nrf2_ich_jloan
        >bash genozip.sh -d "/srv/data/fastq_temp/nrf2_ich_jloan" -p 16

    Here is an example output of the script with only two samples of the nrf2_ich_jloan project:
        >bash genozip.sh -d "/srv/data/fastq_temp/nrf2_ich_jloan"
            Sample folder is set to: /srv/data/fastq_temp/nrf2_ich_jloan
            The following files will be genozipped:
            --------------------------
            /srv/data/fastq_temp/nrf2_ich_jloan/10_Flox_Ma_ICH/10_Flox_Ma_ICH_R1_merged.fastq.gz
            /srv/data/fastq_temp/nrf2_ich_jloan/10_Flox_Mg_ICH/10_Flox_Mg_ICH_R1_merged.fastq.gz
            --------------------------
            Please double check. Do you wish to continue?
            1) Yes
            2) No
            #? 1
            yes
            Compressing with genozip...
            bash -c genozip_fastqgz /srv/data/fastq_temp/nrf2_ich_jloan/10_Flox_Mg_ICH/10_Flox_Mg_ICH_R1_merged.fastq.gz false
            bash -c genozip_fastqgz /srv/data/fastq_temp/nrf2_ich_jloan/10_Flox_Ma_ICH/10_Flox_Ma_ICH_R1_merged.fastq.gz false
            Each sample should have a .fastq and .genozip.fastq file generated.
            The .fastq file is uncompress from the original .fastq.gz file.
            The .genozip.fastq file is uncompress from the .genozip file.
            These two files should have the same md5 hash.
            Checking md5 hash for the following fastq.gz files...
            --------------------------
            /srv/data/fastq_temp/nrf2_ich_jloan/10_Flox_Ma_ICH/10_Flox_Ma_ICH_R1_merged.fastq
            /srv/data/fastq_temp/nrf2_ich_jloan/10_Flox_Ma_ICH/10_Flox_Ma_ICH_R1_merged.genozip.fastq
            /srv/data/fastq_temp/nrf2_ich_jloan/10_Flox_Mg_ICH/10_Flox_Mg_ICH_R1_merged.fastq
            /srv/data/fastq_temp/nrf2_ich_jloan/10_Flox_Mg_ICH/10_Flox_Mg_ICH_R1_merged.genozip.fastq
            --------------------------
            Every md5 has appear exactly twice. okay!
            Deleting the original fastq.gz file:
            --------------------------
            /srv/data/fastq_temp/nrf2_ich_jloan/10_Flox_Ma_ICH/10_Flox_Ma_ICH_R1_merged.fastq.gz
            /srv/data/fastq_temp/nrf2_ich_jloan/10_Flox_Mg_ICH/10_Flox_Mg_ICH_R1_merged.fastq.gz
            --------------------------
            Please check the files to be deleted!!!!! The delete action is **NON-RECOVERABLE**!!! Do you wish to continue?
            1) Yes
            2) No
            #? 1
            yes
            Please **DOUBLE CHECK** the files to be deleted!!!!! The delete action is **NON-RECOVERABLE**!!! Do you wish to continue?
            1) Yes
            2) No
            #? 1
            yes
            Deleting the original fastq.gz file...
            rm /srv/data/fastq_temp/nrf2_ich_jloan/10_Flox_Ma_ICH/10_Flox_Ma_ICH_R1_merged.fastq.gz
            rm /srv/data/fastq_temp/nrf2_ich_jloan/10_Flox_Mg_ICH/10_Flox_Mg_ICH_R1_merged.fastq.gz
            Deleting all fastq (original/genozip) file:
            --------------------------
            /srv/data/fastq_temp/nrf2_ich_jloan/10_Flox_Ma_ICH/10_Flox_Ma_ICH_R1_merged.fastq
            /srv/data/fastq_temp/nrf2_ich_jloan/10_Flox_Ma_ICH/10_Flox_Ma_ICH_R1_merged.genozip.fastq
            /srv/data/fastq_temp/nrf2_ich_jloan/10_Flox_Mg_ICH/10_Flox_Mg_ICH_R1_merged.fastq
            /srv/data/fastq_temp/nrf2_ich_jloan/10_Flox_Mg_ICH/10_Flox_Mg_ICH_R1_merged.genozip.fastq
            --------------------------
            Please check the files to be deleted!!!!! The delete action is **NON-RECOVERABLE**!!! Do you wish to continue?
            1) Yes
            2) No
            #? 1
            yes
            Please **DOUBLE CHECK** the files to be deleted!!!!! The delete action is **NON-RECOVERABLE**!!! Do you wish to continue?
            1) Yes
            2) No
            #? 1
            yes
            Deleting all fastq (original/genozip) file...
            rm /srv/data/fastq_temp/nrf2_ich_jloan/10_Flox_Ma_ICH/10_Flox_Ma_ICH_R1_merged.fastq
            rm /srv/data/fastq_temp/nrf2_ich_jloan/10_Flox_Ma_ICH/10_Flox_Ma_ICH_R1_merged.genozip.fastq
            rm /srv/data/fastq_temp/nrf2_ich_jloan/10_Flox_Mg_ICH/10_Flox_Mg_ICH_R1_merged.fastq
            rm /srv/data/fastq_temp/nrf2_ich_jloan/10_Flox_Mg_ICH/10_Flox_Mg_ICH_R1_merged.genozip.fastq
            Script finished. Size reduced from 2.6G to 1.3G

         >ls -l ./*
            -rw-rw-r-- 1 xinhe xinhe 480 Sep 30 20:36 ./check.md5
            ./10_Flox_Ma_ICH:
            total 668M
            -rw-rw-r-- 1 xinhe xinhe 668M Sep 30 20:35 10_Flox_Ma_ICH_R1_merged.fastq.genozip
            ./10_Flox_Mg_ICH:
            total 609M
            -rw-rw-r-- 1 xinhe xinhe 609M Sep 30 20:35 10_Flox_Mg_ICH_R1_merged.fastq.genozip


    Usage:
      $(basename $0)
            -d folder which contains the sample folders
            -f genozip will overwrite existing genozip file without warning
            -p number of cores. default: 16
            -y disable user input so the script can run programmatically. ** use with caution as script will detele files without asking the user to confirm it.**
            -h see help page
    Example:
        bash /home/xinhe/Projects/sidbdri-utils/scripts/fastqz_to_genozip.sh -d "/srv/data/fastq_temp/nrf2_ich_jloan" -p 16

EOT
}

set -a
function genozip_fastqgz {
    local file=$1
    local force=${2:-false}

    local para=''
    [[ "${force}" == "true" ]] && para='-f'

    para=${para}" --quiet "

    if [  -f "${file}" ]; then
        # extract the fastq.gz to calculate its md5 hash
        gunzip --quiet --keep ${file}
        # genozip the fastq.gz file
        genozip ${para} ${file}
        # after compress with genozip, we uncompress it to check its md5 hash to make sure it is the same as the original one
        genounzip ${para} ${file%.gz}.genozip -o ${file%.fastq.gz}.genozip.fastq
    else
        echo "fastq.gz file <${fastq_gz_file}> does not exist."
        return 1
    fi
}
set +a

## default value
SAMPLE_DIR=''
NUM_CORES=16
CURRENT_DIR=`pwd`
OVERWRITE_EXISTING=false
ENABLE_USER_INPUT=true

## read parameter
while getopts 'd:p:fy' opt; do
  case ${opt} in
    d)  SAMPLE_DIR=${OPTARG}    ;;
    p)  NUM_CORES=${OPTARG}     ;;
    f)  OVERWRITE_EXISTING="true"   ;;
    y)  ENABLE_USER_INPUT="false"   ;;
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



[[ "${SAMPLE_DIR}" == '' ]] && echo "Error: SAMPLE_DIR is not defined." && usage && exit 1
[[ ! "${SAMPLE_DIR}" =~ ^/srv/data/.* ]] && echo "Error: Sampe DIR can only be under /srv/data/" && exit 1
[[ ! -d "${SAMPLE_DIR}"  ]] && echo "Error: SAMPLE_DIR <${SAMPLE_DIR}> does not exist." && usage && exit 1

echo ""
echo "Sample folder is set to: ${SAMPLE_DIR}"
[[ "${OVERWRITE_EXISTING}" == "true" ]] && echo "genozip overwrite is enable. genozip will overwrite existing .genozip."


cd ${SAMPLE_DIR}

echo "The following files will be genozipped:"
echo "--------------------------"
find ${SAMPLE_DIR} -type f -name "*.fastq.gz" | sort
echo "--------------------------"

[[ "${ENABLE_USER_INPUT}" == "true" ]] && confirm

size_before_compress=`du -h -d 1 . | tail -n -1 | cut -f 1`

echo "Compressing with genozip..."
for file in `find ${SAMPLE_DIR} -type f -name "*.fastq.gz"`; do
    echo "genozip_fastqgz ${file} ${OVERWRITE_EXISTING}"
done | xargs -t -n 1 -P ${NUM_CORES} -I % bash -c "%"

echo ""
echo "Each sample should have a .fastq and .genozip.fastq file generated."
echo "The .fastq file is uncompress from the original .fastq.gz file."
echo "The .genozip.fastq file is uncompress from the .genozip file."
echo "These two files should have the same md5 hash."
echo "Checking md5 hash for the following fastq.gz files..."
echo "--------------------------"
find ${SAMPLE_DIR} -type f -name "*.fastq" | sort
echo "--------------------------"


find ${SAMPLE_DIR} -type f -name "*.fastq" | sort | xargs  -P ${NUM_CORES} -n 1 md5sum | sort > check.md5
err=`cat check.md5 | awk '{print $1}' | sort | uniq -c | awk '{if($1!=2){print $2}}'`

if [ "${err}" == '' ]; then
    echo "Every md5 has appear exactly twice. okay! "
else
    echo "The follow md5 hash appear only once, they are different between genozip/original fastq file..."
    echo "${err}"
    exit
fi


echo ""
echo "Deleting the original fastq.gz file:"
echo "--------------------------"
find ${SAMPLE_DIR} -type f -name "*.fastq.gz" | sort
echo "--------------------------"
[[ "${ENABLE_USER_INPUT}" == "true" ]] && confirm "Please check the files to be deleted!!!!! The delete action is **NON-RECOVERABLE**!!! Do you wish to continue?"
[[ "${ENABLE_USER_INPUT}" == "true" ]] && confirm "Please **DOUBLE CHECK** the files to be deleted!!!!! The delete action is **NON-RECOVERABLE**!!! Do you wish to continue?"

echo ""
echo "Deleting the original fastq.gz file..."
find ${SAMPLE_DIR} -type f -name "*.fastq.gz" | sort | xargs -P ${NUM_CORES} -n 1 -t rm

echo ""
echo "Deleting all fastq (original/genozip) file:"
echo "--------------------------"
find ${SAMPLE_DIR} -type f -name "*.fastq" | sort
echo "--------------------------"
[[ "${ENABLE_USER_INPUT}" == "true" ]] && confirm "Please check the files to be deleted!!!!! The delete action is **NON-RECOVERABLE**!!! Do you wish to continue?"
[[ "${ENABLE_USER_INPUT}" == "true" ]] && confirm "Please **DOUBLE CHECK** the files to be deleted!!!!! The delete action is **NON-RECOVERABLE**!!! Do you wish to continue?"

echo ""
echo "Deleting all fastq (original/genozip) file..."
find ${SAMPLE_DIR} -type f -name "*.fastq" | sort | xargs -P ${NUM_CORES} -n 1 -t rm

size_after_compress=`du -h -d 1 . | tail -n -1 | cut -f 1`

echo "Chnaging file group to sidbdri..."
chgrp -Rf sidbdri "${SAMPLE_DIR}"

echo "Script finished. Disk usage from ${size_before_compress} --> ${size_after_compress}"

#find /srv/data/ghardingham/nrf2_ich_jloan -type d | tail -n +2 | xargs -P 31 -n 1 -t -I % bash -c "cp -R % ."