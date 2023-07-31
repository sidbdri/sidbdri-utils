#!/usr/bin/env bash

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
    ["mouse"]="GRCm39"
    ["rat"]="mRatBN7.2"
    ["macaque"]="Mmul_8.0.1"
    ["chimpanzee"]="CHIMP2.1.4"
    ["castaneus"]="CAST_EiJ_v1"
    ["pig"]="Sscrofa11.1"
    ["cow"]="UMD3.1"
)

#check version/url here:
# https://www.ensembl.org/info/website/archives/index.html
declare -A BIOMART_URL=(
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

    wget --user=anonymous --password=${EMAIL} ftp://ftp.ensembl.org/${FILE}
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
