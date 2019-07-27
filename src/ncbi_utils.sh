#!/usr/bin/env bash


function download_srr {
    local out_dir=$1
    local srr=$2
    fastq-dump --split-3 --gzip --outdir ${out_dir}/$srr $srr
}


## @example: download_srrs /home/xinhe/tmp/liver 7 SRR998141 SRR998527 SRR998155 SRR998529 SRR998531 SRR998540 SRR998542
function download_srrs {
    local out_dir=$1
    local nc=$2
    shift 2

    local srrs="$@"

    for srr in ${srrs}; do
        echo "download_srr ${out_dir} ${srr} && fastqc --outdir=${out_dir} ${srr}.fastq.gz 2>${out_dir}/${srr}.log"
    done | xargs -t -n 1 -P ${nc} -I % bash -c "%"
}

function run_fastqc {
    local dir=$1
    local nc=${2:-1}
    for r in `find ${dir} -name '*fastq.gz'`; do
        echo fastqc --outdir=${r%/*} ${r}
    done | xargs -t -n 1 -P ${nc} -I % bash -c "%"
}



