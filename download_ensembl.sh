#!/bin/bash

set -o nounset
set -o errexit

source functions.sh

SPECIES=$1
VERSION=$2
EMAIL=$3

declare -A SCIENTIFIC_NAME=(
    ["mouse"]="mus_musculus"
    ["rat"]="rattus_norvegicus"
)

declare -A ASSEMBLY=(
    ["mouse"]="GRCm38"
    ["rat"]="Rnor_6.0"
)

declare -A ASSEMBLY_TYPE=(
    ["mouse"]="primary_assembly"
    ["rat"]="toplevel"
)

declare -A GENE_DATABASE=(
    ["mouse"]="mmusculus_gene_ensembl"
    ["rat"]="rnorvegicus_gene_ensembl"
)

declare -A BIOMART_URL=(
    ["82"]="www.ensembl.org"
    ["80"]="may2015.archive.ensembl.org"
)

function download_from_ensemble {
    FILE=$1

    wget --user=anonymous --password=${EMAIL} ftp://ftp.ensembl.org/${FILE}
}

scientific_name=${SCIENTIFIC_NAME["$SPECIES"]}
assembly=${ASSEMBLY["$SPECIES"]}
assembly_type=${ASSEMBLY_TYPE["$SPECIES"]}
gene_database=${GENE_DATABASE["$SPECIES"]}
biomart_url=${BIOMART_URL["$VERSION"]}

# Download genome sequence FASTA files

mkdir -p ${assembly_type}
cd ${assembly_type}

genome_fasta=${scientific_name^}.${assembly}.dna.${assembly_type}.fa

download_from_ensemble pub/release-${VERSION}/fasta/${scientific_name}/dna/${genome_fasta}.gz

gunzip ${genome_fasta}.gz
sed -i 's/^>\(.*\) dna.*$/>\1/' ${genome_fasta}
awk '/^>/ {OUT=substr($0,2) ".fa"}; OUT {print >OUT}' ${genome_fasta}
mv ${genome_fasta} ../${SPECIES}_${assembly_type}.fa
cd ..

# Download annotation GTF file

gtf_file=${scientific_name^}.${assembly}.${VERSION}.gtf

download_from_ensemble pub/release-${VERSION}/gtf/${scientific_name}/${gtf_file}.gz

gunzip ${gtf_file}.gz

# Create STAR indices

star_index_dir=STAR_indices/${assembly_type}

mkdir -p ${star_index_dir}

STAR --runThreadN 8 --runMode genomeGenerate --genomeDir ${star_index_dir} --genomeFastaFiles $(list_files ' ' ${assembly_type}/*.fa) --sjdbGTFfile ${gtf_file} --sjdbOverhang 100

# Download gene and ortholog information

wget -O genes.tsv "http://${biomart_url}/biomart/martservice?query=<?xml version=\"1.0\" encoding=\"UTF-8\"?><!DOCTYPE Query><Query  virtualSchemaName = \"default\" formatter = \"TSV\" header = \"0\" uniqueRows = \"0\" count = \"\" datasetConfigVersion = \"0.6\" ><Dataset name = \"${gene_database}\" interface = \"default\" ><Attribute name = \"ensembl_gene_id\" /><Attribute name = \"description\" /><Attribute name = \"chromosome_name\" /><Attribute name = \"external_gene_name\" /></Dataset></Query>"

if [[ "${SPECIES}" != "mouse" ]]; then
    wget -O mouse_orthologs.tsv "http://${biomart_url}/biomart/martservice?query=<?xml version=\"1.0\" encoding=\"UTF-8\"?><!DOCTYPE Query><Query  virtualSchemaName = \"default\" formatter = \"TSV\" header = \"0\" uniqueRows = \"0\" count = \"\" datasetConfigVersion = \"0.6\" ><Dataset name = \"${gene_database}\" interface = \"default\" ><Filter name = \"with_homolog_mmus\" excluded = \"0\"/><Attribute name = \"ensembl_gene_id\" /><Attribute name = \"mmusculus_homolog_ensembl_gene\" /><Attribute name = \"mmusculus_homolog_orthology_type\" /></Dataset></Query>"
fi