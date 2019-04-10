#!/bin/bash
# Download and process data for a particular species from an Ensembl release, e.g.
#
# download_ensembl mouse 95 <your_email> /srv/data/ensembl
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
OUTPUT_DIR=$4

function cleanup {
   echo "Killing all sub-processes..."
   kill -- -$$
}

trap exit INT
trap cleanup EXIT

BASE_DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd -P )
source "${BASE_DIR}/../includes.sh"

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

NUM_THREADS=16

scientific_name=`get_scientific_name ${SPECIES}`
assembly=`get_assembly ${SPECIES}`
assembly_type=`get_assembly_type ${SPECIES}`
biomart_url==`get_biomart_url ${VERSION}`
gene_database=`get_gene_database ${scientific_name}`

# Download genome sequence FASTA files

mkdir -p ${OUTPUT_DIR}
cd ${OUTPUT_DIR}

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

# Create Bowtie indexes

for bowtie2 in ${BOWTIE2_VERSIONS[@]}
do
    bowtie2_version="$(echo ${bowtie2} | sed 's/bowtie2-//')"
    bowtie2_index_dir=BOWTIE2_indices/${assembly_type}_${bowtie2_version}

    mkdir -p ${bowtie2_index_dir}

    bowtie2-build${bowtie2_version} --threads ${NUM_THREADS} *.fa ${bowtie2_index_dir}/bt2index

done

rm -rf BOWTIE2_indices/${assembly_type}
ln -s BOWTIE2_indices/${assembly_type}_${BOWTIE2_VERSIONS[-1]} BOWTIE2_indices/${assembly_type}

# Create Bisulfite index
bismark_genome_preparation --bowtie2 --parallel ${NUM_THREADS} .



# Download gene and ortholog information
download_gene_tb ${SPECIES} ${VERSION} > genes.tsv
download_transcript_tb ${SPECIES} ${VERSION} > transcripts.tsv

if [[ "${SPECIES}" != "mouse" ]]; then
    download_orthologs mouse ${SPECIES} ${VERSION}
fi

if [[ "${SPECIES}" != "human" ]]; then
    download_orthologs human ${SPECIES} ${VERSION}
fi

# Generating refFlat file for Picard RNA-seq metrics
PICARD_DATA=picard
ref_flat=${PICARD_DATA}/${RFF_FILES["$SPECIES"]}
mkdir -p ${PICARD_DATA}

gtfToGenePred -genePredExt -geneNameAsName2 ${gtf_file} ${PICARD_DATA}/refFlat.tmp.txt
paste <(cut -f 12 ${PICARD_DATA}/refFlat.tmp.txt) <(cut -f 1-10 ${PICARD_DATA}/refFlat.tmp.txt) > ${ref_flat}
rm ${PICARD_DATA}/refFlat.tmp.txt
