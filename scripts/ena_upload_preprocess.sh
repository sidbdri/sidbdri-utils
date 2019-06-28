#!/usr/bin/env bash
#
# The script will pre-process sample read files into a format for uploading to
# ENA. The procedure is as follows:
#
# - generate md5sums before copying the read files to output folder
# - copy the read files to the output folder
# - check the md5sums after copying
# - for each sample:
#   - find all read1/read2 files
#   - check read1/read2 have the same number of files
#   - for read1/read2
#     - unzip and concatenate into one file
#     - check the number of lines in the merged file is the sum of pre-merged
#     read files
#     - compress the merged file
#     - move compressed, merged file to output folder
#     - check the compressed, merged file has the same number of lines as the
#     uncompressed, merged file
#   - check merged, compressed read1/read2 files have the same number of rows
# - create a file in the output folder containing md5sums for the merged,
# compressed read files
# - finish

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
The script will pre-process sample read files into a format for uploading to
ENA. The procedure is as follows:

- generate md5sums before copying the read files to output folder
- copy the read files to the output folder
- check the md5sums after copying
- for each sample:
  - find all read1/read2 files
  - check read1/read2 have the same number of files
  - for read1/read2
    - unzip and concatenate into one file
    - check the number of lines in the merged file is the sum of pre-merged
    read files
    - compress the merged file
    - move compressed, merged file to output folder
    - check the compressed, merged file has the same number of lines as the
    uncompressed, merged file
  - check merged, compressed read1/read2 files have the same number of rows
- create a file in the output folder containing md5sums for the merged,
compressed read files
- finish

Usage:
  $(basename $0)
        -i input folder which contains the sample folders
        -o output folder which will contain the read files for upload
        -s space separated sample names to be included, e.g. "1N1 1N2 1N3" 
        (default: all samples in input folder)
        -1 pattern for read 1 files. (default: "*_1.fastq.gz")
        -2 pattern for read 1 files. (default: "*_2.fastq.gz")
        -p number of cores. (default: 1)
        -t test-mode, which will use only the first 10 lines from read files.
        (default: false)
        -h display help
EOT
}

## default values

SAMPLE_NAMES=''
READ1_IDENTIFIER="*_1.fastq.gz"
READ2_IDENTIFIER="*_2.fastq.gz"
NUM_CORES=1
TEST_MODE=false

while getopts 'i:o:s:1:2:p:t' opt; do
  case ${opt} in
    i)  DATA_DIR=${OPTARG}    ;;
    o)  OUT_DIR=${OPTARG}    ;;
    s)  SAMPLE_NAMES=${OPTARG}    ;;
    1)  READ1_IDENTIFIER=${OPTARG}        ;;
    2)  READ2_IDENTIFIER=${OPTARG}            ;;
    p)  NUM_CORES=${OPTARG}            ;;
    t)  TEST_MODE=true      ;;
    *)  usage; exit         ;;
  esac
done

if [ -z "${SAMPLE_NAMES}" ]; then
    echo "-s is not set, using all sampleis found in ${DATA_DIR}."
    SAMPLE_NAMES=`ls ${DATA_DIR}| grep -v .txt | paste -sd " " -`
    echo ${SAMPLE_NAMES}
fi

FINAL_DIR=${OUT_DIR}/final

# check data folder exists and is not empty
if [  ! -d "${DATA_DIR}" ]; then
    >&2 echo "Error: data folder <${DATA_DIR}> does not exist." && exit 1
fi

# check all samples exist
for sample in ${SAMPLE_NAMES}; do
    if [  ! -d "${DATA_DIR}/${sample}" ]; then
        >&2 echo "Error: sample folder <${DATA_DIR}/${sample}> does not exist." && exit 1
    fi
done

# check output folder does not exist
if [  -d "${OUT_DIR}" ]; then
    >&2 echo "Error: output folder <${OUT_DIR}> already exists." && exit 1
fi
mkdir -p ${OUT_DIR} ${FINAL_DIR}

# generate md5sums before copying the read files to output folder
echo "Generate md5sums before copying the read files to output folder..."
mkdir -p ${OUT_DIR}/md5sum
for sample in ${SAMPLE_NAMES}; do
    echo "cd ${DATA_DIR} && find ${sample} -type f -exec md5sum {} + > ${OUT_DIR}/md5sum/${sample}.chk"
done | xargs -t -n 1 -P ${NUM_CORES} -I % bash -c "%"

# copy the samples to output folder
echo "Copying reads to output folder..."
for sample in ${SAMPLE_NAMES}; do
    echo "cp -R -t ${OUT_DIR} ${DATA_DIR}/${sample}"
done | xargs -t -n 1 -P ${NUM_CORES} -I % bash -c "%"

# check md5sums after copying
echo "Checking md5sums of read files in output folder..."
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
    # If any invocation of the command exits with a status of 255,
    # xargs will stop immediately without reading any further input.
    # An error message is issued on stderr when this happens
    tmp_dir=$1
    sample=$2
    read_index=$3
    out_dir=$4
    test_mode=$5
    read_files=${@:6}

    merged_file=${tmp_dir}/${sample}/merged/${sample}_${read_index}.fastq
    merged_compress_file=${merged_file}.gz
    out_file=${out_dir}/${sample}_${read_index}.fastq.gz
    rm -rf ${merged_file} ${merged_compress_file} ${out_file}

    # unzip the reads and concat into one file
    if [ "$test_mode" = true ]; then
        for f in ${read_files}; do
            zcat -cdk ${f} | head -10 >>  ${merged_file}
        done
    else
        for f in ${read_files}; do
           gzip -cdk ${f} >>  ${merged_file}
         done
    fi

    # check if merged reads file has the same number of lines as the sum of
    # the pre-merged read files
    lines_in_merged_file=`wc -l ${merged_file} | awk '{print $1}'`

    if [ "$test_mode" = true ]; then
       lines_in_pre_merged_file=`for f in ${read_files}; do
                                 zcat ${f} | head -10 |wc -l
                              done | awk '{ SUM += $0} END { print SUM }'`
    else
        lines_in_pre_merged_file=`for f in ${read_files}; do
                                 zcat ${f} | wc -l
                              done | awk '{ SUM += $0} END { print SUM }'`
    fi

    if [ ! ${lines_in_merged_file} -eq ${lines_in_pre_merged_file} ];then
         >&2 echo "Error: number of rows differ after concatenating read files.  \
         Expecting: ${lines_in_pre_merged_file}    " `basename ${merged_file}`": ${lines_in_merged_file} "
         exit 255
    fi

    # compress the merged read files
    gzip -fk ${merged_file}

    # move the compressed file to the output folder
    mv ${merged_compress_file} ${out_file}

    # check the compressed merged read files has the same number of rows as the
    # uncompressed file
    lines_in_merged_compressed_out_file=`zcat ${out_file} | wc -l`
    if [ ! ${lines_in_merged_file} -eq ${lines_in_merged_compressed_out_file} ];then
         >&2 echo "Error: number of rows differ after compressing. \
         Expecting : ${lines_in_merged_file}     " `basename ${merged_file}`": ${lines_in_merged_compressed_out_file} "
         exit 255
    fi
}

function check_two_pairs_equal_lines(){
    ## check two fastq.gz files have the same number of rows
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

    echo "merge_reads ${OUT_DIR} ${sample} 1 ${FINAL_DIR} ${TEST_MODE} " ${read_1_files}
    echo "merge_reads ${OUT_DIR} ${sample} 2 ${FINAL_DIR} ${TEST_MODE} " ${read_2_files}
done |  xargs -t -n 1 -P ${NUM_CORES} -I % bash -c "%"

# check if the two final fastq.gz files for each sample have the same number of
# rows
echo "Checking final reads has the same number of rows between pairs..."
for sample in ${SAMPLE_NAMES}; do
    echo "check_two_pairs_equal_lines ${FINAL_DIR}/${sample}_1.fastq.gz ${FINAL_DIR}/${sample}_2.fastq.gz"
done | xargs -t -n 1 -P ${NUM_CORES} -I % bash -c "%"

# generate md5sums for the final fastq.gz files
echo "Generate md5sum for final reads..."
for sample in ${SAMPLE_NAMES}; do
    echo "cd ${FINAL_DIR} && md5sum ${sample}_1.fastq.gz"
    echo "cd ${FINAL_DIR} && md5sum ${sample}_2.fastq.gz"
done | xargs -t -n 1 -P ${NUM_CORES} -I % bash -c "%" > ${FINAL_DIR}/md5sum.txt

echo "Finished! The files for upload can be found in ${FINAL_DIR}"

