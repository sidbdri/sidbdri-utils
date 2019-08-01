#!/usr/bin/env bash

# Returns files matching the specified pattern(s) separated by a particular
# delimiter character.
function list_files {
    local DELIMITER=$1
    shift
    local FILES=$@

    LIST=$(ls -1 $FILES | tr '\n' "$DELIMITER")
    echo ${LIST%$DELIMITER}
}

function sam2bam {
    local sam=$1
    local bam=$2

    sambamba view -S ${sam} -f bam > ${bam}
}

function subset_fastq_gz {
    local file=$1
    local number_lines=$2
    local out=$3
    zcat ${file} | head -${2} | gzip > ${out}
}

function pf {
    local file=$1
    local number_of_lines=${2:-10}



    if [[ ${file} == *.bam ]]; then
        sambamba view ${file} | head -${number_of_lines}
    fi

    if [[ ${file} == *.sam ]]; then
        cat ${file} | awk ' /^[^@]/ {print $0}' | head -${number_of_lines}
    fi

}

function grep_from_bam {
    local file=$1
    local pattern=$2

    sambamba view ${file} | grep "${pattern}"
}

#subset_fastq_with_reads_from_bam SRR5467502_1.fastq.gz 10 /home/xinhe/tmp/SRR5467502_test/filtered_reads/SRR5467502___mouse___filtered.bam /home/xinhe/tmp/SRR5467502_test/filtered_reads/SRR5467502___rat___filtered.bam
function subset_fastq_with_reads_from_bam {
    local fastq=$1
    local number_line_from_bam=$2
    local bam=$3

    pf ${bam} ${number_line_from_bam} | cut -f 1 | sort -u | sed 's/^/@/g'  | sed 's/$/ /g'
    zcat ${fastq} | grep -A 3 `pf ${bam} ${number_line_from_bam} | cut -f 1 | sort -u | sed 's/^/@/g'  | sed 's/$/ /g'`
}

