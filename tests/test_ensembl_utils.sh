#!/usr/bin/env bash

source "$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd -P )/../includes.sh"

set -o nounset
set -o errexit

get_assembly 'human' '103'

get_assembly_type 'huamn'

get_gtf_file 'human' '95'

get_gene_database 'homo_sapiens'

get_scientific_name 'human'

get_biomart_url 95

get_primary_chromosomes human 95 | head -1

download_orthologs human mouse 95 | head -1

query_biomart 95 '<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE Query>
<Query  virtualSchemaName = "default" formatter = "TSV" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" >

	<Dataset name = "hsapiens_gene_ensembl" interface = "default" >
		<Attribute name = "ensembl_gene_id" />
		<Attribute name = "ensembl_transcript_id" />
	</Dataset>
</Query>' true true | head -2

download_gene_tb 'human' 95 | head -1

download_transcript_tb 'human' 95 | head -1
