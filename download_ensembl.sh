#!/bin/bash

# Download and process data for a particular species from an Ensembl release, e.g.
#
# download_ensembl mouse 82 <your_email>
#
# The script:
# 1) Downloads top-level sequences for the species' genome in FASTA format
# (primary assembly sequences for human and mouse),
# 2) Downloads gene annotations in GTF format,
# 3) Creates a genome index for the STAR mapping tool,
# 4) Builds a Salmon index for the transcripts defined in the GTF file
# 5) Downloads gene information (gene id, name, chromosome, description) into a
# file 'genes.tsv'. Only genes on primary assembly sequences are retained (i.e.
# those on patch sequences are removed),
# 6) Download information on mouse and human orthologous genes to files
# '<species>_orthologs.tsv'.

set -o nounset
set -o errexit

source functions.sh

SPECIES=$1
VERSION=$2
EMAIL=$3

declare -A SCIENTIFIC_NAME=(
    ["human"]="homo_sapiens"
    ["mouse"]="mus_musculus"
    ["rat"]="rattus_norvegicus"
)

declare -A ASSEMBLY=(
    ["human"]="GRCh38"
    ["mouse"]="GRCm38"
    ["rat"]="Rnor_6.0"
)

declare -A BIOMART_URL=(
    ["85"]="www.ensembl.org"
    ["84"]="mar2016.archive.ensembl.org"
    ["83"]="dec2015.archive.ensembl.org"
    ["82"]="sep2015.archive.ensembl.org"
    ["81"]="jul2015.archive.ensembl.org"
    ["80"]="may2015.archive.ensembl.org"
)

function download_from_ensembl {
    local FILE=$1

    wget --user=anonymous --password=${EMAIL} ftp://ftp.ensembl.org/${FILE}
}

function get_assembly_type {
    local SPECIES=$1 

    if [ "${SPECIES}" == "human" ] || [ "${SPECIES}" == "mouse" ] ; then
        echo "primary_assembly"
    else
        echo "toplevel"
    fi
}

function get_gene_database {
    local SCIENTIFIC_NAME=$1

    echo $(echo ${scientific_name} | sed 's/\(.\).*_\(.*\)/\1\2/')_gene_ensembl
}

function download_orthologs {
    local ORTHOLOG_SPECIES=$1
    local BIOMART_URL=$2
    local GENE_DATABASE=$3

    ortho_sci_name=${SCIENTIFIC_NAME["$ORTHOLOG_SPECIES"]}
    ortho_short_name=$(echo ${ortho_sci_name} | sed 's/\(.\).*_\(.*\)/\1\2/')
    ortho_shorter_name=${ortho_short_name:0:4}

    wget -O ${ORTHOLOG_SPECIES}_orthologs.tsv "http://${BIOMART_URL}/biomart/martservice?query=<?xml version=\"1.0\" encoding=\"UTF-8\"?><!DOCTYPE Query><Query  virtualSchemaName = \"default\" formatter = \"TSV\" header = \"0\" uniqueRows = \"0\" count = \"\" datasetConfigVersion = \"0.6\" ><Dataset name = \"${GENE_DATABASE}\" interface = \"default\" ><Filter name = \"with_homolog_${ortho_shorter_name}\" excluded = \"0\"/><Attribute name = \"ensembl_gene_id\" /><Attribute name = \"${ortho_short_name}_homolog_ensembl_gene\" /><Attribute name = \"${ortho_short_name}_homolog_orthology_type\" /></Dataset></Query>"
}

NUM_THREADS=16

scientific_name=${SCIENTIFIC_NAME["$SPECIES"]}
assembly=${ASSEMBLY["$SPECIES"]}
assembly_type=$(get_assembly_type ${SPECIES})
biomart_url=${BIOMART_URL["$VERSION"]}
gene_database=$(get_gene_database ${scientific_name})

# Download genome sequence FASTA files

mkdir -p ${assembly_type}
cd ${assembly_type}

genome_fasta=${scientific_name^}.${assembly}.dna.${assembly_type}.fa

download_from_ensembl pub/release-${VERSION}/fasta/${scientific_name}/dna/${genome_fasta}.gz

gunzip ${genome_fasta}.gz
sed -i 's/^>\(.*\) dna.*$/>\1/' ${genome_fasta}
awk '/^>/ {OUT=substr($0,2) ".fa"}; OUT {print >OUT}' ${genome_fasta}
mv ${genome_fasta} ../${SPECIES}_${assembly_type}.fa
cd ..

# Download annotation GTF file
#
gtf_file=${scientific_name^}.${assembly}.${VERSION}.gtf

download_from_ensembl pub/release-${VERSION}/gtf/${scientific_name}/${gtf_file}.gz

gunzip ${gtf_file}.gz

# Create STAR indices

star_index_dir=STAR_indices/${assembly_type}

mkdir -p ${star_index_dir}

STAR --runThreadN ${NUM_THREADS} --runMode genomeGenerate --genomeDir ${star_index_dir} --genomeFastaFiles $(list_files ' ' ${assembly_type}/*.fa) --sjdbGTFfile ${gtf_file} --sjdbOverhang 100

# Create Salmon and Kallisto indexes

salmon_index_dir=salmon_index
kallisto_index_dir=kallisto_index

transcripts_ref=${salmon_index_dir}/transcripts
transcripts_fasta=${transcripts_ref}.transcripts.fa

mkdir -p ${salmon_index_dir}

rsem-prepare-reference -p ${NUM_THREADS} --gtf ${gtf_file} ${assembly_type} ${transcripts_ref}

salmon index -p ${NUM_THREADS} -t ${transcripts_fasta} -i ${salmon_index_dir}
echo "Salmon index created with $(salmon --version 2>&1)." >> README

kallisto index -i ${kallisto_index} ${transcripts_fasta}
echo "Kallisto index created with $(kallisto version)." >> README

# Download gene and ortholog information

wget -qO- "http://${biomart_url}/biomart/martservice?query=<?xml version=\"1.0\" encoding=\"UTF-8\"?><!DOCTYPE Query><Query  virtualSchemaName = \"default\" formatter = \"TSV\" header = \"0\" uniqueRows = \"0\" count = \"\" datasetConfigVersion = \"0.6\" ><Dataset name = \"${gene_database}\" interface = \"default\" ><Attribute name = \"ensembl_gene_id\" /><Attribute name = \"description\" /><Attribute name = \"chromosome_name\" /><Attribute name = \"external_gene_name\" /></Dataset></Query>" |\
    awk -F'\t' 'NR==FNR {a[$0]=$0} NR>FNR {if($3==a[$3]) print $0}' <(ls -1 ${assembly_type} | sed 's/.fa//') - > genes.tsv

if [[ "${SPECIES}" != "mouse" ]]; then
    download_orthologs mouse ${biomart_url} ${gene_database}
fi

if [[ "${SPECIES}" != "human" ]]; then
    download_orthologs human ${biomart_url} ${gene_database}
fi
