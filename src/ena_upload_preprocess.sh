#!/usr/bin/env bash
#
# The script will pre-process the sample read files into
# a format for uploading to ENA. The procedure is as follows:
#
# - generate md5sum before copy the reads to out folder
# - copy the reads to out folder
# - check the md5sum after copy
# - for each sample:
#   - find all read1/read2
#   - check read1/read2 has the same number of files
#   - for read1/read2
#     - unzip and concat into one file
#     - check number of rows in merge file is the sum of pre-merged read files
#     - compress merged file
#     - move to output folder
#     - check the compressed merged file has the same number of rows as the uncompress merged file
#   - check merged-compressed read1/read2 has the same number of rows
# - Finish
#
#
# Note from Owen:
# Multi-file case: /srv/data/ghardingham/dreadd_phasel
# Single-file case: /srv/data/tspiresjones/app_tau
#
# Inputs:
# - directory containing the read files: e.g. /srv/data/ghardingham/dreadd_phasel
# - set of sample names: e.g. Gq_Cnt_A, Gs_4h_CNO_A, Gs_Cnt_A
# - patterns that identify _1 and _2 files
#
# Process:
# - make a copy of original files and work on those, rather than original files
# - for each sample
# - -  for “_1” “_2” files (bearing in mind how the two types are named)
# - - - calculate number of lines in the original files
# - - - unzip and concatenate original files
# - - - count lines in merged file and check equal
# - - - gzip new file and call <sample_name>_1/2.fastq.gz
# - - check that final _1 and _2 files have same number of lines
# - - check that final _1 and _2 files are different
#
# Output
# - a folder containing, for each sample, <sample_name>_1.fastq.gz, <sample_name>_2.fastq.gz
# - text file containing md5sums for the *.fastq.gz

set -o nounset
set -o errexit
#set -o xtrace

function cleanup {
#   echo "Killing all sub-processes..."
   kill -- -$$
}

trap exit INT
trap cleanup EXIT


function usage {
  cat <<EOT
The script will pre-process the sample read files into
a format for uploading to ENA. The procedure is as follows:

- generate md5sum before copy the reads to out folder
- copy the reads to out folder
- check the md5sum after copy
- for each sample:
  - find all read1/read2
  - check read1/read2 has the same number of files
  - for read1/read2
    - unzip and concat into one file
    - check number of rows in merge file is the sum of pre-merged read files
    - compress merged file
    - move to output folder
    - check the compressed merged file has the same number of rows as the uncompress merged file
  - check merged-compressed read1/read2 has the same number of rows
- Finish

Usage:
  $(basename $0)
        -i input folder which contains the sample folders
        -o output folder which contains the read-to-upload reads
        -s space separated sample_names to be included. default: all samples in input folder.
        -1 pattern for read 1 files. default: "*_1.fastq.gz"
        -2 pattern for read 1 files. default: ""*_2.fastq.gz""
        -p number of cores. default: 1
        -h see help page
Example:
    bash ena_upload_preprocess.sh  -i /srv/data/ghardingham/dreadd_phasel -o ~/tmp/dreadd_phasel -s "Gs_Cnt_A Gs_Cnt_B" -1 "*_1.fastq.gz" -2 "*_2.fastq.gz" -p 10

EOT
}


#DATA_DIR=/srv/data/ghardingham/dreadd_phasel

#bash ena_upload_preprocess.sh  -i /srv/data/ghardingham/dreadd_phasel -o ~/tmp/dreadd_phasel_test -s "Gs_Cnt_A" -1 "*_1.fastq.gz" -2 "*_2.fastq.gz" -p 10
#OUT_DIR=~/tmp/${DATA_DIR##*/}

## default value
SAMPLE_NAMES=''
READ1_IDENTIFIER="*_1.fastq.gz"
READ2_IDENTIFIER="*_2.fastq.gz"
NUM_CORES=1

while getopts 'i:o:s:1:2:p:' opt; do
  case ${opt} in
    i)  DATA_DIR=${OPTARG}    ;;
    o)  OUT_DIR=${OPTARG}    ;;
    s)  SAMPLE_NAMES=${OPTARG}    ;;
    1)  READ1_IDENTIFIER=${OPTARG}        ;;
    2)  READ2_IDENTIFIER=${OPTARG}            ;;
    p)  NUM_CORES=${OPTARG}            ;;
    *)  usage; exit         ;;
  esac
done

if [ -z ${SAMPLE_NAMES} ]; then
    echo "-s is not set, using all sample found in ${DATA_DIR}."
    SAMPLE_NAMES=`ls ${DATA_DIR}| grep -v .txt | paste -sd " " -`
    echo ${SAMPLE_NAMES}
fi

#
FINAL_DIR=${OUT_DIR}/final
#SAMPLE_NAMES="Gs_Cnt_A Gs_Cnt_B"


# check data folder exists and not empty
if [  ! -d "${DATA_DIR}" ]; then
    >&2 echo "Error: data folder <${DATA_DIR}> does not exist." && exit 1
fi

# check all samples exist
for sample in ${SAMPLE_NAMES}; do
    if [  ! -d "${DATA_DIR}/${sample}" ]; then
        >&2 echo "Error: sample folder <${DATA_DIR}/${sample}> does not exist." && exit 1
    fi
done

# check out folder does not exist
if [  -d "${OUT_DIR}" ]; then
    >&2 echo "Error: out folder <${OUT_DIR}> already exist." && exit 1
fi
mkdir -p ${OUT_DIR} ${FINAL_DIR}


# generate md5sum before copy the reads to out folder
echo "Generate md5sum before copy the reads to out folder..."
mkdir -p ${OUT_DIR}/md5sum
for sample in ${SAMPLE_NAMES}; do
    echo "cd ${DATA_DIR} && find ${sample} -type f -exec md5sum {} + > ${OUT_DIR}/md5sum/${sample}.chk"
done | xargs -t -n 1 -P ${NUM_CORES} -I % bash -c "%"

# we copy the samples to out folder
echo "Copying reads to out folder..."
for sample in ${SAMPLE_NAMES}; do
    echo "cp -R -t ${OUT_DIR} ${DATA_DIR}/${sample}"
done | xargs -t -n 1 -P ${NUM_CORES} -I % bash -c "%"

# we check the md5sum after copy
echo "Checking md5sum of reads in out folder..."
for sample in ${SAMPLE_NAMES}; do
    echo "cd ${OUT_DIR} && md5sum -c ${OUT_DIR}/md5sum/${sample}.chk"
done | xargs -t -n 1 -P ${NUM_CORES} -I % bash -c "%" > ${OUT_DIR}/md5sum/check1.txt

if grep -q 'FAILED' "${OUT_DIR}/md5sum/check1.txt"; then
    >&2 echo "Error: md5sum check failed. The following file did not pass the check:"
    >&2 grep 'FAILED' "${OUT_DIR}/md5sum/check1.txt"
    exit 1
fi



set -a
function merge_reads(){
    #If any invocation of the command exits with a status of 255,
    # xargs will stop immediately without reading any further input.
    # An error message is issued on stderr when this happens
    tmp_dir=$1
    sample=$2
    read_index=$3
    out_dir=$4
    read_files=${@:5}

    merged_file=${tmp_dir}/${sample}/merged/${sample}_${read_index}.fastq
    merged_compress_file=${merged_file}.gz
    out_file=${out_dir}/${sample}_${read_index}.fastq.gz
    rm -rf ${merged_file} ${merged_compress_file} ${out_file}

    ## we unzip the reads and concat into one file
    for f in ${read_files}; do
       gzip -cdk ${f} >>  ${merged_file}
    done

#    for f in ${read_files}; do
#       zcat -cdk ${f} | head -4 >>  ${merged_file}
#    done

    ## check if merged reads file has the same number of rows as the sum of the pre-merged read files
    lines_in_merged_file=`wc -l ${merged_file} | awk '{print $1}'`
    lines_in_pre_merged_file=`for f in ${read_files}; do
                                 zcat ${f} | wc -l
                              done | awk '{ SUM += $0} END { print SUM }'`

    if [ ! ${lines_in_merged_file} -eq ${lines_in_pre_merged_file} ];then
         >&2 echo "Error: number of rows differ after concatenating read files.  \
         Excepting: ${lines_in_pre_merged_file}    " `basename ${merged_file}`": ${lines_in_merged_file} "
         exit 255
    fi

    ## we compress the merged read files
    gzip -fk ${merged_file}

    ## move the compressed file to the output folder
    mv ${merged_compress_file} ${out_file}

    ## we check the compressed merged read files has the same number of rows as the umcompressed file
    lines_in_merged_compressed_out_file=`zcat ${out_file} | wc -l`
    if [ ! ${lines_in_merged_file} -eq ${lines_in_merged_compressed_out_file} ];then
         >&2 echo "Error: number of rows differ after compressing. \
         Excepting : ${lines_in_merged_file}     " `basename ${merged_file}`": ${lines_in_merged_compressed_out_file} "
         exit 255
    fi

}

function checkTwoPairEqualRow(){
    ## check two fastq.gz has the same number of rows
    read_1=$1
    read_2=$2

    number_of_lines_read_1=`zcat ${read_1} | wc -l`
    number_of_lines_read_2=`zcat ${read_2} | wc -l`

    if [ ! ${number_of_lines_read_1} -eq ${number_of_lines_read_2} ];then
        >&2 echo "Error: number of rows differ between two pairs of reads."\
        ${read_1}" : "${number_of_lines_read_1}"     "${read_2}" : "${number_of_lines_read_2}
        exit 255
    fi
}
set +a


echo "Merging samples..."
for sample in ${SAMPLE_NAMES}; do
    mkdir -p ${OUT_DIR}/${sample}/merged
    read_1_files="$(ls ${OUT_DIR}/${sample}/${READ1_IDENTIFIER})"
    read_2_files="$(ls ${OUT_DIR}/${sample}/${READ2_IDENTIFIER})"

    ## check number_of_read1_files = number_of_read2_files
    number_of_read1_files=`echo ${read_1_files} | tr ' ' '\n' | wc -l`
    number_of_read2_files=`echo ${read_2_files} | tr ' ' '\n' | wc -l`
    if [ ! ${number_of_read1_files} -eq ${number_of_read2_files} ]; then
         >&2 echo "Error: The number of files for read_1 and read_2 is not equal for sample: ${sample}"
         >&2 echo "Read 1 found:"
         >&2 echo "${read_1_files}"
         >&2 echo "Read 2 found:"
         >&2 echo "${read_2_files}"
         exit 1
    fi

    echo "merge_reads ${OUT_DIR} ${sample} 1 ${FINAL_DIR} " ${read_1_files}
    echo "merge_reads ${OUT_DIR} ${sample} 2 ${FINAL_DIR} " ${read_2_files}
done |  xargs -t -n 1 -P ${NUM_CORES} -I % bash -c "%"


## we check if the two final fastq.gz for each sample have the same number of rows
echo "Checking final reads has the same number of rows between pairs..."
for sample in ${SAMPLE_NAMES}; do
    echo "checkTwoPairEqualRow ${FINAL_DIR}/${sample}_1.fastq.gz ${FINAL_DIR}/${sample}_2.fastq.gz"
done | xargs -t -n 1 -P ${NUM_CORES} -I % bash -c "%"

## we generate md5sum for final fastq.gz
echo "Generate md5sum for final reads..."
for sample in ${SAMPLE_NAMES}; do
    echo "md5sum ${FINAL_DIR}/${sample}_1.fastq.gz"
    echo "md5sum ${FINAL_DIR}/${sample}_2.fastq.gz"
done | xargs -t -n 1 -P ${NUM_CORES} -I % bash -c "%" > ${FINAL_DIR}/md5sum.txt

echo "Finished! The ready-to-upload reads can be found ${FINAL_DIR}"

