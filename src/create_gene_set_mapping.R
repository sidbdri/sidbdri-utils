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
msigdb_path <- str_c('/srv/data/genome/', SPECIES, '/msigdb/v2025.1')
out_dir <- file.path(OUT_DIR,'msigdb','v2025.1')
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

file_list <- list(
  "human" = list(
    "CURATED" = "c2.all.v2025.1.Hs.entrez.gmt",
    "MOTIF" = "c3.all.v2025.1.Hs.entrez.gmt",
    "GO" = "c5.all.v2025.1.Hs.entrez.gmt",
    "MSIGDB_CELL_TYPE" = "c8.all.v2025.1.Hs.entrez.gmt"
  ),
  "mouse" = list(
    "CURATED" = "m2.all.v2025.1.Mm.entrez.gmt",
    "MOTIF" = "m3.all.v2025.1.Mm.entrez.gmt",
    "GO" = "m5.all.v2025.1.Mm.entrez.gmt",
    "MSIGDB_CELL_TYPE" = "m8.all.v2025.1.Mm.entrez.gmt"
  )
)

# the referece species is the species annotation file to be use.
# for example, for rat, we are using mouse annotation file, so reference_species is 'mouse'
# for mouse/human, reference_species is the same as species so can leave it as NA
get_gene_sets <- function(species, gene_set_name, reference_species=NA) {
  if (species == 'rat') {
    reference_species <- 'mouse'
  }
  if (is.na(reference_species)) {
    reference_species <- species
  }
  message(str_c('Processing ', species, ' ', gene_set_name, ' using ', reference_species, ' as reference'))
  
  file_name <- file_list[[reference_species]][[gene_set_name]]
  msigdb_data <- file.path(msigdb_path,file_name) %>% read_lines()
  
  gene_set_names <- map_chr(msigdb_data, function(x) {str_split(x, "\t")[[1]][1]})
  gene_sets <- map(msigdb_data, function(x) {str_split(str_split(x, "\t", n=3)[[1]][3], "\t")[[1]]})
  
  names(gene_sets) <- gene_set_names
  
  if (species %in% names(file_list)) {
    ret <- gene_sets
  }else{
    pb <- txtProgressBar(max = length(gene_sets), style = 3)
    count <- 0
    
    ortholog_info <- get_reference_vs_species_ortholog_info(species, reference_species)
    reference_entrez_col <- str_c(reference_species, "_entrez_id")
    
    ret <- gene_sets %>% map(function(y) {
      count <<- count + 1
      setTxtProgressBar(pb, count)
      ortholog_info[which(ortholog_info[[reference_entrez_col]] %in% y),]$species_entrez_id %>% unique
    })
    
    close(pb)
  }
  
  
  saveRDS(ret,file = file.path(out_dir,str_c(gene_set_name,'.all.v2025.1.entrez.gmt.Rdata')) )
}


get_reference_vs_species_ortholog_info <- function(species, reference_species = "human") {
  ortholog_file <- str_c(reference_species, "_orthologs.tsv")
  reference_gene_col <- str_c(reference_species, "_gene")
  reference_entrez_col <- str_c(reference_species, "_entrez_id")
  
  orthologs <- species %>%
    str_c("/srv/data/genome/", ., "/ensembl-", ENSEMBL_VERSION, "/", ortholog_file) %>%
    read_tsv(col_names = c("species_gene", reference_gene_col, "type"))
  
  species_genes <- get_gene_info(species)
  reference_genes <- get_gene_info(reference_species)
  
  orthologs_entrez <- orthologs %>%
    left_join(species_genes %>%
                dplyr::select(gene, entrez_id) %>%
                dplyr::rename(species_gene = gene, species_entrez_id = entrez_id)) %>%
    left_join(reference_genes %>%
                dplyr::select(gene, entrez_id) %>%
                dplyr::rename(!!sym(reference_gene_col) := gene, !!sym(reference_entrez_col) := entrez_id)) %>%
    dplyr::select(species_entrez_id, !!sym(reference_entrez_col)) %>%
    filter(!is.na(species_entrez_id) & !is.na(!!sym(reference_entrez_col))) %>%
    distinct()
  
  same_name_entrez <- (species_genes %>%
                         mutate(gene_name = tolower(gene_name)) %>%
                         dplyr::select(gene_name, entrez_id) %>%
                         filter(!is.na(entrez_id) & !is.na(gene_name)) %>%
                         dplyr::rename(species_entrez_id = entrez_id)) %>%
    inner_join(reference_genes %>%
                 mutate(gene_name = tolower(gene_name)) %>%
                 dplyr::select(gene_name, entrez_id) %>%
                 filter(!is.na(entrez_id) & !is.na(gene_name)) %>%
                 dplyr::rename(!!sym(reference_entrez_col) := entrez_id)) %>%
    dplyr::select(-gene_name) %>%
    distinct()
  
  orthologs_entrez %>%
    rbind(same_name_entrez) %>%
    distinct()
}


get_gene_info <- function(species) {
  species %>%
    str_c("/srv/data/genome/", ., "/ensembl-",ENSEMBL_VERSION,"/genes.tsv") %>%
    read_tsv(col_names = c("gene", "description", "chromosome", "gene_name", "entrez_id", "gene_type"),
             col_types = list(chromosome = col_character())) %>%
    group_by(gene) %>%
    dplyr::filter(row_number()==1) %>%
    ungroup
}


gene_set_categories <- list("CURATED", "MOTIF", "GO", "MSIGDB_CELL_TYPE")

gene_set_categories %>%
  set_names(.) %>%
  lapply(function(category, ...) get_gene_sets(SPECIES, category)) %>% invisible()