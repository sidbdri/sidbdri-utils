# Rscript xx.R mouse 100 /srv/data/genome/mouse/ensembl-100
library(dplyr)
library(magrittr)
library(readr)
library(stringr)
library(purrr)

args = commandArgs(trailingOnly=TRUE)

if (length(args)!=3) {
    stop("Example: Rscript human 100 /srv/data/genome/mouse/ensembl-100")
}

SPECIES <- args[1]
ENSEMBL_VERSION <- args[2]
OUT_DIR <- args[3]

root<-str_c('/srv/data/genome/', SPECIES, '/ensembl-',ENSEMBL_VERSION,'/')
msigdb_path <- '/srv/data/genome/human/msigdb/v7.0'
out_dir <- file.path(OUT_DIR,'msigdb','v7.0')
if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
}


get_gene_sets <- function(species, gene_set_name) {
    msigdb_data <- file.path(msigdb_path,str_c(gene_set_name,'.all.v7.0.entrez.gmt')) %>%
    read_lines()

    gene_set_names <- map_chr(msigdb_data, function(x) {str_split(x, "\t")[[1]][1]})
    gene_sets <- map(msigdb_data, function(x) {str_split(str_split(x, "\t", n=3)[[1]][3], "\t")[[1]]})

    names(gene_sets) <- gene_set_names

    if (species == "human") {
        ret <- gene_sets
    }else{
        pb <- txtProgressBar(max = length(gene_sets), style = 3)
        count <- 0

        ortholog_info <- get_human_vs_species_ortholog_info(species)

        ret <- gene_sets %>% map(function(y) {
            count <<- count + 1
            setTxtProgressBar(pb, count)
            ortholog_info[which(ortholog_info$human_entrez_id %in% y),]$species_entrez_id %>% unique
        })

        close(pb)
    }


    saveRDS(ret,file = file.path(out_dir,str_c(gene_set_name,'.all.v7.0.entrez.gmt.Rdata')) )
}


get_human_vs_species_ortholog_info <- function(species) {
    orthologs <- species %>%
        str_c("/srv/data/genome/", ., "/ensembl-",ENSEMBL_VERSION,"/human_orthologs.tsv") %>%
        read_tsv(col_names = c("species_gene", "human_gene", "type"))

    species_genes <- get_gene_info(species)
    human_genes <- get_gene_info("human")

    orthologs_entrez <- orthologs %>%
        left_join(species_genes %>%
            dplyr::select(gene, entrez_id) %>%
            dplyr::rename(species_gene = gene, species_entrez_id = entrez_id)) %>%
        left_join(human_genes %>%
            dplyr::select(gene, entrez_id) %>%
            dplyr::rename(human_gene = gene, human_entrez_id = entrez_id)) %>%
        dplyr::select(species_entrez_id, human_entrez_id) %>%
        filter(!is.na(species_entrez_id) & !is.na(human_entrez_id)) %>%
        distinct()

    same_name_entrez <- (species_genes %>%
        mutate(gene_name = tolower(gene_name)) %>%
        dplyr::select(gene_name, entrez_id) %>%
        filter(!is.na(entrez_id) & !is.na(gene_name)) %>%
        dplyr::rename(species_entrez_id = entrez_id)) %>%
        inner_join(human_genes %>%
            mutate(gene_name = tolower(gene_name)) %>%
            dplyr::select(gene_name, entrez_id) %>%
            filter(!is.na(entrez_id) & !is.na(gene_name)) %>%
            dplyr::rename(human_entrez_id = entrez_id)) %>%
        dplyr::select(-gene_name) %>%
        distinct()

    human_species_entrez_mappings <- orthologs_entrez %>%
        rbind(same_name_entrez) %>%
        distinct
}


get_gene_info <- function(species) {
    species %>%
        str_c("/srv/data/genome/", ., "/ensembl-",ENSEMBL_VERSION,"/genes.tsv") %>%
        read_tsv(col_names = c("gene", "description", "chromosome", "gene_name", "entrez_id", "gene_type"),
        col_types = list(chromosome = col_character())) %>%
        group_by(gene) %>%
        filter(row_number()==1) %>%
        ungroup
}





gene_set_categories <- list("CURATED", "MOTIF", "GO", "CELL_TYPE")

gene_set_categories %>%
    set_names(.) %>%
    lapply(function(category, ...) get_gene_sets(SPECIES, category)) %>% invisible()



