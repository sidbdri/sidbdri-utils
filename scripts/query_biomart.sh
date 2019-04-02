#!/usr/bin/env bash

function usage {
  cat <<EOT

Usage:
  $(basename $0) -s human [-v 94]
        -v ensembl version. defaul: latest
        -q query xml file.
        -p peek the first 10 lines of the query results
        -h see help page
Example:
    bash query_biomart.sh  -s human -v latest -p -q ./gene_transcript_APPRIS.xml

EOT
}


VERSION="95"
SPECIES='human'
PEEK=false


while getopts 'v:q:p' opt; do
  case ${opt} in
    v)  VERSION=${OPTARG}    ;;
    q)  XML=${OPTARG}        ;;
    p)  PEEK=true            ;;
    *)  usage; exit         ;;
  esac
done

BASE_DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd -P )

source "${BASE_DIR}/../includes.sh"

query_biomart ${VERSION} "${XML}" true ${PEEK}







