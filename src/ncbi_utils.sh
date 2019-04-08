#!/usr/bin/env bash


function download_srr {
    local out_dir=$1
    local srr=$2
    fastq-dump --split-3 --gzip --outdir ${out_dir} $srr
}


function download_srrs {
    local out_dir=$1; shift
    local nc=$2; shift
    local srrs=$@

    for srr in ${srrs}; do
        echo "download_srr ${out_dir} ${srr}"
    done| xargs -t -n 1 -P ${nc} -I % bash -c "%"
}
