#!/usr/bin/env bash

declare -A SCIENTIFIC_NAME=(
    ["human"]="homo_sapiens"
    ["mouse"]="mus_musculus"
    ["mouse_cba"]="mus_musculus_cbaj"
    ["rat"]="rattus_norvegicus"
    ["macaque"]="macaca_mulatta"
    ["chimpanzee"]="pan_troglodytes"
    ["castaneus"]="mus_musculus_casteij"
    ["pig"]="sus_scrofa"
    ["cow"]="bos_taurus"
)

declare -A ASSEMBLY=(
    ["human"]="GRCh38"
    ["mouse"]="GRCm39"
    ["mouse_cba"]="CBA_J_v1"
    ["rat"]="GRCr8"
    ["macaque"]="Mmul_8.0.1"
    ["chimpanzee"]="CHIMP2.1.4"
    ["castaneus"]="CAST_EiJ_v1"
    ["pig"]="Sscrofa11.1"
    ["cow"]="UMD3.1"
)

#check version/url here:
# https://www.ensembl.org/info/website/archives/index.html
declare -A BIOMART_URL=(
    ["116"]="jun2026.archive.ensembl.org"
    ["115"]="sep2025.archive.ensembl.org"
#    ["114"]="may2025.archive.ensembl.org"
    ["113"]="oct2024.archive.ensembl.org"
    ["112"]="may2024.archive.ensembl.org"
    ["111"]="jan2024.archive.ensembl.org"
    ["110"]="jul2023.archive.ensembl.org"
    ["109"]="feb2023.archive.ensembl.org"
    ["108"]="oct2022.archive.ensembl.org"
    ["107"]="jul2022.archive.ensembl.org"
    ["106"]="apr2022.archive.ensembl.org"
    ["105"]="dec2021.archive.ensembl.org"
    ["104"]="may2021.archive.ensembl.org"
    ["103"]="feb2021.archive.ensembl.org"
    ["102"]="nov2020.archive.ensembl.org"
    ["101"]="aug2020.archive.ensembl.org"
    ["100"]="apr2020.archive.ensembl.org"
    ["99"]="jan2020.archive.ensembl.org/"
    ["98"]="sep2019.archive.ensembl.org/"
    ["97"]="jul2019.archive.ensembl.org/"
    ["96"]="apr2019.archive.ensembl.org"
    ["95"]="jan2019.archive.ensembl.org"
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
    ["84"]="mar2016.archive.ensembl.org"
    ["83"]="dec2015.archive.ensembl.org"
    ["82"]="sep2015.archive.ensembl.org"
    ["81"]="jul2015.archive.ensembl.org"
    ["80"]="may2015.archive.ensembl.org"
)


function download_from_ensembl {
    local FILE=$1

    wget --user=anonymous --password=${EMAIL} http://ftp.ensembl.org/${FILE}
}

function get_assembly {
    local SPECIES=$1
    local VERSION=$2

    if [ "${SPECIES}" == "mouse" ] && [ "${VERSION}" -le "102" ] ; then
        echo "GRCm38"
    elif [ "${SPECIES}" == "rat" ] && [ "${VERSION}" -le "104" ] ; then
        echo "Rnor_6.0"
    else
        echo ${ASSEMBLY["$SPECIES"]}
    fi
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

    local scientific_name=`get_scientific_name ${SPECIES}`
    local assembly=`get_assembly ${SPECIES} ${VERSION}`

    if [ "${SPECIES}" == "castaneus" ] ; then
        VERSION=86
    fi

    echo ${scientific_name^}.${assembly}.${VERSION}.gtf
}

function get_rff_file {
    local SPECIES=$1
    local VERSION=$2

    local scientific_name=`get_scientific_name ${SPECIES}`
    local assembly=`get_assembly ${SPECIES} ${VERSION}`

    echo ${scientific_name^}.${assembly}.${VERSION}.rff
}

function get_gene_database {
    local SCIENTIFIC_NAME=$1

    echo $(echo ${SCIENTIFIC_NAME} | sed 's/\(.\).*_\(.*\)/\1\2/')_gene_ensembl
}

function get_scientific_name {
    local SPECIES=$1

    echo ${SCIENTIFIC_NAME["$SPECIES"]}
}

function get_biomart_url {
    local VERSION=$1

    echo ${BIOMART_URL["${VERSION}"]}
}

function get_primary_chromosomes {
    local SPECIES=$1
    local VERSION=$2

    local enmsebl_folder="/srv/data/genome"
    local assembly_type=`get_assembly_type ${SPECIES}`

    ls -1 ${enmsebl_folder}/${SPECIES}/ensembl-${VERSION}/${assembly_type} | sed 's/.fa//'
}

function download_orthologs {
    local SPECIES=$1
    local ORTHOLOG_SPECIES=$2
    local ENSEMBL_VERSION=$3

    local scientific_name=`get_scientific_name ${SPECIES}`
    local gene_database=`get_gene_database ${scientific_name}`
    local assembly_type=`get_assembly_type ${SPECIES}`

    local ortho_sci_name=`get_scientific_name ${ORTHOLOG_SPECIES}`
    local ortho_short_name=$(echo ${ortho_sci_name} | sed 's/\(.\).*_\(.*\)/\1\2/')
    local ortho_shorter_name=${ortho_short_name:0:4}

    if [ "${ENSEMBL_VERSION}" -ge "86" ]; then
        filter_name="with_${ortho_short_name}_homolog"
    else
        filter_name="with_homolog_${ortho_shorter_name}"
    fi

    local query=`cat <<EOT
<?xml version="1.0" encoding="UTF-8"?>
    <!DOCTYPE Query>
        <Query  virtualSchemaName = "default" formatter = "TSV" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" >
            <Dataset name = "${gene_database}" interface = "default" >
            <Filter name = "${filter_name}" excluded = "0"/>
            <Attribute name = "ensembl_gene_id" />
            <Attribute name = "${ortho_short_name}_homolog_ensembl_gene" />
            <Attribute name = "${ortho_short_name}_homolog_orthology_type" />
    </Dataset>
</Query>
EOT`

    query_biomart ${ENSEMBL_VERSION} "${query}" false false
}

function query_biomart {
    local VERSION=$1
    local XML=$2
    local HEADER=${3:-true}
    local PEEK=${4:-false}

    if [  -s "${XML}" ]; then
        local query="$(cat ${XML})"
        local header=`grep 'Attribute name' ${XML} | sed -E 's/.+= \"(.+)\".+/\1/' | tr "\n" "\t" | sed 's/\t$/\n/'`
    else
        local query="${XML}"
        local header=`echo ${query} | sed -E 's/>/>\n/g' | grep 'Attribute name' | sed -E 's/.+= \"(.+)\".+/\1/' | tr "\n" "\t" | sed 's/\t$/\n/'`
    fi

    local biomart_url=`get_biomart_url ${VERSION}`

    query="http://${biomart_url}/biomart/martservice?query="$(echo $query | tr -d '\n')
    >&2 echo "query: "$query

    if [ "${HEADER}" = true ] ; then
        echo "${header}"
    fi

    if [ "${PEEK}" = true ] ; then
        wget -qO- "$query" | head
    else
        wget -qO- "$query"
    fi
}

function download_gene_tb {
    local SPECIES=$1
    local VERSION=$2

    local scientific_name=`get_scientific_name ${SPECIES}`
    local gene_database=`get_gene_database ${scientific_name}`
    local assembly_type=`get_assembly_type ${SPECIES}`

    #in version 96, entrezgene attribute is changed to entrezgene_id
    entrezgene_str=entrezgene
    if [ ${VERSION} -ge 97 ];then
        entrezgene_str=entrezgene_id
    fi

    query=`echo '<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE Query>
<Query  virtualSchemaName = "default" formatter = "TSV" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" >
	<Dataset name = "gene_database" interface = "default" >
		<Attribute name = "ensembl_gene_id" />
		<Attribute name = "description" />
		<Attribute name = "chromosome_name" />
		<Attribute name = "external_gene_name" />
		<Attribute name = "entrezgene_str" />
		<Attribute name = "gene_biotype" />
</Dataset>
</Query>
' |  sed "s/gene_database/${gene_database}/" | sed "s/entrezgene_str/${entrezgene_str}/"`

#     filter out non-primary chromosomes
    query_biomart ${VERSION} "${query}" false false | \
    awk -F'\t' 'NR==FNR {a[$0]=$0} NR>FNR {if($3==a[$3]) print $0}' <( get_primary_chromosomes ${SPECIES} ${VERSION} ) -
}

function download_transcript_tb {
    local SPECIES=$1
    local VERSION=$2

    local scientific_name=`get_scientific_name ${SPECIES}`
    local gene_database=`get_gene_database ${scientific_name}`
    local assembly_type=`get_assembly_type ${SPECIES}`

    query=`echo '<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE Query>
	<Query  virtualSchemaName = "default" formatter = "TSV" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" >
		<Dataset name = "gene_database" interface = "default" >
			<Attribute name = "ensembl_transcript_id" />
			<Attribute name = "transcript_biotype" />
			<Attribute name = "ensembl_gene_id" />
			<Attribute name = "chromosome_name" />
</Dataset>
</Query>
' |  sed "s/gene_database/${gene_database}/"`

    # filter out non-primary chromosomes
    query_biomart ${VERSION} "${query}" false false | \
    awk -F'\t' 'NR==FNR {a[$0]=$0} NR>FNR {if($4==a[$4]) print $0}' <( get_primary_chromosomes ${SPECIES} ${VERSION} ) -
}


function generate_picard_refFlat  {
    local PICARD_DATA_DIR=$1
    local SPECIES=$2
    local VERSION=$3
    local GTF_FILE=$4

    local ref_flat=`get_rff_file ${SPECIES} ${VERSION}`

    mkdir -p ${PICARD_DATA_DIR}

    gtfToGenePred -genePredExt -geneNameAsName2 ${GTF_FILE} ${PICARD_DATA_DIR}/refFlat.tmp.txt
    paste <(cut -f 12 ${PICARD_DATA_DIR}/refFlat.tmp.txt) <(cut -f 1-10 ${PICARD_DATA_DIR}/refFlat.tmp.txt) > ${PICARD_DATA_DIR}/${ref_flat}
    rm ${PICARD_DATA_DIR}/refFlat.tmp.txt
}

# Verify that download_ensembl produced the expected genome/annotation/index files.
# Uses STAR_VERSIONS / SALMON_VERSIONS / BOWTIE2_VERSIONS from the caller when set.
# Usage: verify_ensembl_download <species> <version> <output_dir>
function verify_ensembl_download {
    local SPECIES=$1
    local VERSION=$2
    local OUTPUT_DIR=$3

    local assembly_type=`get_assembly_type ${SPECIES}`
    local gtf_file=$(get_gtf_file ${SPECIES} ${VERSION})
    local rff_file=$(get_rff_file ${SPECIES} ${VERSION})
    local genome_fasta=${SPECIES}_${assembly_type}.fa

    local star_versions=("${STAR_VERSIONS[@]}")
    local salmon_versions=("${SALMON_VERSIONS[@]}")
    local bowtie2_versions=("${BOWTIE2_VERSIONS[@]}")
    if [ ${#star_versions[@]} -eq 0 ]; then star_versions=("STAR2.7.11b"); fi
    if [ ${#salmon_versions[@]} -eq 0 ]; then salmon_versions=("salmon1.10.0"); fi
    if [ ${#bowtie2_versions[@]} -eq 0 ]; then bowtie2_versions=("bowtie2-2.3.4.3" "bowtie2-2.4.4"); fi

    local missing=()
    local failed=0

    _check_exists() {
        local path=$1
        local desc=${2:-$1}
        if [ ! -e "${OUTPUT_DIR}/${path}" ]; then
            missing+=("MISSING: ${desc} (${path})")
            failed=1
            return 1
        fi
        return 0
    }

    _check_nonempty() {
        local path=$1
        local desc=${2:-$1}
        if [ ! -s "${OUTPUT_DIR}/${path}" ]; then
            missing+=("EMPTY/MISSING: ${desc} (${path})")
            failed=1
            return 1
        fi
        return 0
    }

    _check_min_size() {
        local path=$1
        local min_bytes=$2
        local desc=${3:-$1}
        local full="${OUTPUT_DIR}/${path}"
        if [ ! -s "${full}" ]; then
            missing+=("EMPTY/MISSING: ${desc} (${path})")
            failed=1
            return 1
        fi
        local size
        size=$(stat -c%s "${full}" 2>/dev/null || stat -f%z "${full}")
        if [ "${size}" -lt "${min_bytes}" ]; then
            missing+=("TOO SMALL (${size} < ${min_bytes} bytes): ${desc} (${path})")
            failed=1
            return 1
        fi
        return 0
    }

    # Biomart occasionally returns an HTML error page; reject that.
    _check_tsv() {
        local path=$1
        local desc=${2:-$1}
        local full="${OUTPUT_DIR}/${path}"
        if [ ! -s "${full}" ]; then
            missing+=("EMPTY/MISSING: ${desc} (${path})")
            failed=1
            return 1
        fi
        if head -c 200 "${full}" | grep -qiE '<(!DOCTYPE|html)|Service unavailable'; then
            missing+=("NOT A TSV (looks like HTML/error page): ${desc} (${path})")
            failed=1
            return 1
        fi
        if ! grep -q $'\t' "${full}"; then
            missing+=("NO TAB SEPARATORS: ${desc} (${path})")
            failed=1
            return 1
        fi
        return 0
    }

    # tx2gene is space-separated (awk default), not tab-separated
    _check_two_col() {
        local path=$1
        local desc=${2:-$1}
        local full="${OUTPUT_DIR}/${path}"
        if [ ! -s "${full}" ]; then
            missing+=("EMPTY/MISSING: ${desc} (${path})")
            failed=1
            return 1
        fi
        if ! awk 'NF!=2 {bad=1; exit} END{exit bad+0}' "${full}"; then
            missing+=("EXPECTED 2 COLUMNS: ${desc} (${path})")
            failed=1
            return 1
        fi
        return 0
    }

    _check_symlink() {
        local path=$1
        local desc=${2:-$1}
        local full="${OUTPUT_DIR}/${path}"
        if [ ! -L "${full}" ]; then
            missing+=("MISSING SYMLINK: ${desc} (${path})")
            failed=1
            return 1
        fi
        if [ ! -e "${full}" ]; then
            missing+=("BROKEN SYMLINK: ${desc} (${path} -> $(readlink "${full}"))")
            failed=1
            return 1
        fi
        return 0
    }

    echo "=============================================="
    echo "Verifying Ensembl download: ${SPECIES} ${VERSION}"
    echo "Output dir: ${OUTPUT_DIR}"
    echo "=============================================="

    # Genome FASTA + per-chromosome split
    _check_min_size "${genome_fasta}" 100000000 "genome FASTA"
    _check_exists "${assembly_type}" "assembly chromosome dir"
    local n_chr
    n_chr=$(find "${OUTPUT_DIR}/${assembly_type}" -maxdepth 1 -name '*.fa' 2>/dev/null | wc -l | tr -d ' ')
    if [ "${n_chr}" -lt 20 ]; then
        missing+=("TOO FEW chromosome FASTAs in ${assembly_type}/ (${n_chr} < 20)")
        failed=1
    fi

    # GTF
    _check_min_size "${gtf_file}" 10000000 "GTF annotation"

    # STAR indices
    local star star_version star_index_dir
    for star in "${star_versions[@]}"; do
        star_version="$(echo ${star} | sed 's/STAR//')"
        star_index_dir=STAR_indices/${assembly_type}_${star_version}
        _check_min_size "${star_index_dir}/Genome" 100000000 "STAR Genome (${star_version})"
        _check_min_size "${star_index_dir}/SA" 100000000 "STAR SA (${star_version})"
        _check_nonempty "${star_index_dir}/SAindex" "STAR SAindex (${star_version})"
        _check_nonempty "${star_index_dir}/genomeParameters.txt" "STAR genomeParameters (${star_version})"
    done
    _check_symlink "STAR_indices/${assembly_type}" "STAR default symlink"

    # RSEM transcript reference + Salmon indices
    _check_min_size "transcripts_ref/transcripts.transcripts.fa" 10000000 "RSEM transcripts FASTA"
    _check_nonempty "transcripts_ref/transcripts.ti" "RSEM transcripts.ti"
    local salmon salmon_version salmon_index_dir
    for salmon in "${salmon_versions[@]}"; do
        salmon_version="$(echo ${salmon} | sed 's/salmon//')"
        salmon_index_dir=SALMON_indices/${assembly_type}_${salmon_version}
        _check_nonempty "${salmon_index_dir}/info.json" "Salmon info.json (${salmon_version})"
        _check_min_size "${salmon_index_dir}/seq.bin" 1000000 "Salmon seq.bin (${salmon_version})"
    done
    _check_symlink "SALMON_indices/${assembly_type}" "Salmon default symlink"

    # Bowtie2 indices
    local bowtie2 bowtie2_version bowtie2_index_dir bt2_suffix
    for bowtie2 in "${bowtie2_versions[@]}"; do
        bowtie2_version="$(echo ${bowtie2} | sed 's/bowtie2-//')"
        bowtie2_index_dir=BOWTIE2_indices/${assembly_type}_${bowtie2_version}
        for bt2_suffix in 1.bt2 2.bt2 3.bt2 4.bt2 rev.1.bt2 rev.2.bt2; do
            _check_nonempty "${bowtie2_index_dir}/bt2index.${bt2_suffix}" "Bowtie2 ${bowtie2_version} bt2index.${bt2_suffix}"
        done
    done
    _check_symlink "BOWTIE2_indices/${assembly_type}" "Bowtie2 default symlink"

    # Bismark Bisulfite genome
    for bt2_suffix in 1.bt2 2.bt2 3.bt2 4.bt2 rev.1.bt2 rev.2.bt2; do
        _check_nonempty "Bisulfite_Genome/CT_conversion/BS_CT.${bt2_suffix}" "Bisulfite CT BS_CT.${bt2_suffix}"
        _check_nonempty "Bisulfite_Genome/GA_conversion/BS_GA.${bt2_suffix}" "Bisulfite GA BS_GA.${bt2_suffix}"
    done

    # Gene / transcript tables
    _check_tsv "genes.tsv" "genes.tsv"
    _check_tsv "transcripts.tsv" "transcripts.tsv"
    _check_two_col "tx2gene.tsv" "tx2gene.tsv"
    _check_nonempty "gene_lengths.csv" "gene_lengths.csv"

    if [[ "${SPECIES}" != "mouse" ]]; then
        _check_tsv "mouse_orthologs.tsv" "mouse_orthologs.tsv"
    fi
    if [[ "${SPECIES}" != "human" ]]; then
        _check_tsv "human_orthologs.tsv" "human_orthologs.tsv"
    fi

    # Picard refFlat
    _check_min_size "picard/${rff_file}" 1000000 "Picard refFlat"

    # Gene-set mappings produced by create_gene_set_mapping.R
    local gs
    for gs in CURATED MOTIF GO MSIGDB_CELL_TYPE; do
        _check_nonempty "msigdb/v2025.1/${gs}.all.v2025.1.entrez.gmt.Rdata" "msigdb ${gs}"
    done

    echo
    if [ "${failed}" -eq 0 ]; then
        echo "SUCCESS: all expected files present for ${SPECIES} ensembl-${VERSION}"
        return 0
    fi

    echo "FAILURE: download/index verification failed for ${SPECIES} ensembl-${VERSION}"
    printf '  %s\n' "${missing[@]}"
    return 1
}
