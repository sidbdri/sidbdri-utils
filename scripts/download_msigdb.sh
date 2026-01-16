#!/usr/bin/env bash
set -e

VERSION="2025.1"

download_msigdb() {
    local species=$1
    local species_code=$2
    local dest="/srv/data/genome/${species}/msigdb/v${VERSION}"
    
    mkdir -p "$dest"
    cd "$dest"
    
    local url="https://data.broadinstitute.org/gsea-msigdb/msigdb/release/${VERSION}.${species_code}/msigdb_v${VERSION}.${species_code}_files_to_download_locally.zip"
    local zipfile="msigdb_v${VERSION}.${species_code}_files_to_download_locally.zip"
    
    curl -L -O "$url"
    unzip -o -j "$zipfile"
}

download_msigdb "human" "Hs"
download_msigdb "mouse" "Mm"
