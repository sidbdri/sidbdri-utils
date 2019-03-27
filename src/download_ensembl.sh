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

SPECIES=$1
VERSION=$2
EMAIL=$3

function cleanup {
   echo "Killing all sub-processes..."
   kill -- -$$
}

trap exit INT
trap cleanup EXIT

declare -A SCIENTIFIC_NAME=(
    ["human"]="homo_sapiens"
    ["mouse"]="mus_musculus"
    ["rat"]="rattus_norvegicus"
    ["macaque"]="macaca_mulatta"
    ["chimpanzee"]="pan_troglodytes"
    ["castaneus"]="mus_musculus_casteij"
    ["pig"]="sus_scrofa"
    ["cow"]="bos_taurus"
)

declare -A ASSEMBLY=(
    ["human"]="GRCh38"
    ["mouse"]="GRCm38"
    ["rat"]="Rnor_6.0"
    ["macaque"]="Mmul_8.0.1"
    ["chimpanzee"]="CHIMP2.1.4"
    ["castaneus"]="CAST_EiJ_v1"
    ["pig"]="Sscrofa11.1"
    ["cow"]="UMD3.1"
)

declare -A BIOMART_URL=(
    ["95"]="www.ensembl.org"
    ["94"]="oct2018.archive.ensembl.org"
    ["93"]="jul2018.archive.ensembl.org"
    ["92"]="apr2018.archive.ensembl.org"
    ["91"]="dec2017.archive.ensembl.org"
    ["90"]="aug2017.archive.ensembl.org"
    ["89"]="may2017.archive.ensembl.org"
    ["88"]="mar2017.archive.ensembl.org"
    ["87"]="dec2016.archive.ensembl.org"
    ["86"]="oct2016.archive.ensembl.org"
    ["85"]="jul2016.archive.ensembl.org"
    ["84"]="mar2016.archive.ensembl.org"REF_FLAT
    ["83"]="dec2015.archive.ensembl.org"
    ["82"]="sep2015.archive.ensembl.org"
    ["81"]="jul2015.archive.ensembl.org"
    ["80"]="may2015.archive.ensembl.org"
)

STAR_VERSIONS=(
    "STAR2.5.3a"
    "STAR2.6.1d"
)

KALLISTO_VERSIONS=(
    "kallisto0.43.1"
    "kallisto0.45.0"
)

RSEM_VERSION=1.3.1

SALMON_VERSIONS=(
    "salmon0.8.2"
    "salmon0.12.0"
)

BOWTIE2_VERSIONS=(
    "bowtie2-2.3.4"
    "bowtie2-2.3.4.3"
)

declare -A RFF_FILES=(
        ["mouse"]="Mus_musculus.GRCm38.${VERSION}.rff"
        ["human"]="Homo_sapiens.GRCh38.${VERSION}.rff"
        ["rat"]="Rattus_norvegicus.Rnor_6.0.${VERSION}.rff"
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

function get_gtf_file {
    local SPECIES=$1
    local VERSION=$2

    scientific_name=${SCIENTIFIC_NAME["$SPECIES"]}
    assembly=${ASSEMBLY["$SPECIES"]}

    if [ "${SPECIES}" == "castaneus" ] ; then
        VERSION=86
    fi

    echo ${scientific_name^}.${assembly}.${VERSION}.gtf
}

function get_gene_database {
    local SCIENTIFIC_NAME=$1

    echo $(echo ${scientific_name} | sed 's/\(.\).*_\(.*\)/\1\2/')_gene_ensembl
}

function download_orthologs {
    local ORTHOLOG_SPECIES=$1
    local BIOMART_URL=$2
    local GENE_DATABASE=$3
    local ENSEMBL_VERSION=$4

    ortho_sci_name=${SCIENTIFIC_NAME["$ORTHOLOG_SPECIES"]}
    ortho_short_name=$(echo ${ortho_sci_name} | sed 's/\(.\).*_\(.*\)/\1\2/')
    ortho_shorter_name=${ortho_short_name:0:4}

    if [ "${ENSEMBL_VERSION}" -ge "86" ]; then
        filter_name="with_${ortho_short_name}_homolog"
    else
        filter_name="with_homolog_${ortho_shorter_name}"
    fi

    wget -O ${ORTHOLOG_SPECIES}_orthologs.tsv "http://${BIOMART_URL}/biomart/martservice?query=<?xml version=\"1.0\" encoding=\"UTF-8\"?> <!DOCTYPE Query> <Query  virtualSchemaName = \"default\" formatter = \"TSV\" header = \"0\" uniqueRows = \"0\" count = \"\" datasetConfigVersion = \"0.6\" > <Dataset name = \"${GENE_DATABASE}\" interface = \"default\" > <Filter name = \"${filter_name}\" excluded = \"0\"/> <Attribute name = \"ensembl_gene_id\" /> <Attribute name = \"${ortho_short_name}_homolog_ensembl_gene\" /> <Attribute name = \"${ortho_short_name}_homolog_orthology_type\" /> </Dataset> </Query>"
}

# Returns files matching the specified pattern(s) separated by a particular
# delimiter character.
function list_files {
    local DELIMITER=$1
    shift
    local FILES=$@
        
    LIST=$(ls -1 $FILES | tr '\n' "$DELIMITER")
    echo ${LIST%$DELIMITER}    
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

gtf_file=$(get_gtf_file ${SPECIES} ${VERSION})

download_from_ensembl pub/release-${VERSION}/gtf/${scientific_name}/${gtf_file}.gz  

gunzip -c ${gtf_file}.gz | tail -n +6 > ${gtf_file}
rm ${gtf_file}.gz



# Create STAR indices
for star in ${STAR_VERSIONS[@]}
do
    star_version="$(echo ${star} | sed 's/STAR//')"
    star_index_dir=${assembly_type}_${star_version}

    mkdir -p STAR_indices/${star_index_dir}

    ${star} --runThreadN ${NUM_THREADS} --runMode genomeGenerate --genomeDir STAR_indices/${star_index_dir} --genomeFastaFiles $(list_files ' ' ${assembly_type}/*.fa) --sjdbGTFfile ${gtf_file} --sjdbOverhang 100
done
rm -rf STAR_indices/${assembly_type}
ln -s STAR_indices/${assembly_type}_${STAR_VERSIONS[-1]} STAR_indices/${assembly_type}


# Create Salmon and Kallisto indexes
transcripts_ref=transcripts_ref/transcripts
transcripts_fasta=${transcripts_ref}.transcripts.fa
mkdir -p ${transcripts_ref}
rsem-prepare-reference${RSEM_VERSION} -p ${NUM_THREADS} --gtf ${gtf_file} ${assembly_type} ${transcripts_ref}

for salmon in ${SALMON_VERSIONS[@]}
do
    salmon_version="$(echo ${salmon} | sed 's/salmon//')"
    salmon_index_dir=SALMON_indices/${assembly_type}_${salmon_version}

    mkdir -p ${salmon_index_dir}

    salmon index -p ${NUM_THREADS} -t ${transcripts_fasta} -i ${salmon_index_dir}
    echo "Salmon index created with $(salmon --version 2>&1)." >> README
done
rm -rf SALMON_indices/${assembly_type}
ln -s SALMON_indices/${assembly_type}_${SALMON_VERSIONS[-1]} SALMON_indices/${assembly_type}

for kallisto in ${KALLISTO_VERSIONS[@]}
do
    kallisto_version="$(echo ${kallisto} | sed 's/kallisto//')"
    kallisto_index_dir=KALLISTO_indices/${assembly_type}_${kallisto_version}

    mkdir -p ${kallisto_index_dir}

    kallisto index -i ${kallisto_index_dir}/kallisto_index ${transcripts_fasta}
    echo "Kallisto index created with $(kallisto --version 2>&1)." >> README
done
rm -rf KALLISTO_indices/${assembly_type}
ln -s KALLISTO_indices/${assembly_type}_${KALLISTO_VERSIONS[-1]}/kallisto_index KALLISTO_indices/${assembly_type}


# Create bowtie indexes

for bowtie2 in ${BOWTIE2_VERSIONS[@]}
do
    bowtie2_version="$(echo ${bowtie2} | sed 's/bowtie2-//')"
    bowtie2_index_dir=BOWTIE2_indices/${assembly_type}_${bowtie2_version}

    mkdir -p ${bowtie2_index_dir}

    bowtie2-build${bowtie2_version} --threads ${NUM_THREADS} *.fa ${bowtie2_index_dir}/bt2index

done
rm -rf BOWTIE2_indices/${assembly_type}
ln -s BOWTIE2_indices/${assembly_type}_${BOWTIE2_VERSIONS[-1]} BOWTIE2_indices/${assembly_type}

# Download gene and ortholog information

wget -qO- "http://${biomart_url}/biomart/martservice?query=<?xml version=\"1.0\" encoding=\"UTF-8\"?><!DOCTYPE Query><Query  virtualSchemaName = \"default\" formatter = \"TSV\" header = \"0\" uniqueRows = \"0\" count = \"\" datasetConfigVersion = \"0.6\" ><Dataset name = \"${gene_database}\" interface = \"default\" ><Attribute name = \"ensembl_gene_id\" /><Attribute name = \"description\" /><Attribute name = \"chromosome_name\" /><Attribute name = \"external_gene_name\" /><Attribute name = \"entrezgene\" /><Attribute name = \"gene_biotype\" /></Dataset></Query>" |\
    awk -F'\t' 'NR==FNR {a[$0]=$0} NR>FNR {if($3==a[$3]) print $0}' <(ls -1 ${assembly_type} | sed 's/.fa//') - > genes.tsv

wget -qO- "http://${biomart_url}/biomart/martservice?query=<?xml version=\"1.0\" encoding=\"UTF-8\"?><!DOCTYPE Query><Query  virtualSchemaName = \"default\" formatter = \"TSV\" header = \"0\" uniqueRows = \"0\" count = \"\" datasetConfigVersion = \"0.6\" ><Dataset name = \"${gene_database}\" interface = \"default\" ><Attribute name = \"ensembl_transcript_id\" /><Attribute name = \"transcript_biotype\" /><Attribute name = \"ensembl_gene_id\" /><Attribute name = \"chromosome_name\" /></Dataset></Query>" |\
    awk -F'\t' 'NR==FNR {a[$0]=$0} NR>FNR {if($4==a[$4]) print $0}' <(ls -1 ${assembly_type} | sed 's/.fa//') - > transcripts.tsv

if [[ "${SPECIES}" != "mouse" ]]; then
    download_orthologs mouse ${biomart_url} ${gene_database} ${VERSION}
fi

if [[ "${SPECIES}" != "human" ]]; then
    download_orthologs human ${biomart_url} ${gene_database} ${VERSION}
fi

# Generating refFlat file for Picard RNA-seq metrics
PICARD_DATA=picard
ref_flat=${PICARD_DATA}/${RFF_FILES["$SPECIES"]}
mkdir -p ${PICARD_DATA}
gtfToGenePred -genePredExt -geneNameAsName2 ${gtf_file} ${PICARD_DATA}/refFlat.tmp.txt
paste <(cut -f 12 ${PICARD_DATA}/refFlat.tmp.txt) <(cut -f 1-10 ${PICARD_DATA}/refFlat.tmp.txt) > ${ref_flat}
rm ${PICARD_DATA}/refFlat.tmp.txt



