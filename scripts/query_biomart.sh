#!/usr/bin/env bash

VERSION="95"
PEEK=false

function usage {
  cat <<EOT

Usage:
  $(basename $0)
        -v ensembl version. default: ${VERSION}
        -q query xml file.
        -p peek the first 10 lines of the query results
        -h see help page
Example:
    bash query_biomart.sh -v 94 -p -q ./gene_transcript_APPRIS.xml

    bash query_biomart.sh -v 94 -p -q '<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE Query>
<Query  virtualSchemaName = "default" formatter = "TSV" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" >
    <Dataset name = "hsapiens_gene_ensembl" interface = "default" >
        <Attribute name = "ensembl_gene_id" />
        <Attribute name = "external_gene_name" />
    </Dataset>
</Query>'

EOT
}

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







