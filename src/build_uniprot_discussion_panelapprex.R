# # build_uniprot_discussion_panelapprex.R
# 
# library(dplyr)
# library(tidyr)
# library(readr)
# library(rtracklayer)
# library(UniprotR)
# 
# output_dir <- "../data"
# 
# dir.create(file.path(output_dir, "ontology_taxa"), showWarnings = FALSE, recursive = TRUE)
# dir.create(file.path(output_dir, "uniprot_panelapp"), showWarnings = FALSE, recursive = TRUE)
# 
# # paths for UniProt source files
# path_uniprot_gff <- "../ref/uniprot/uniprot-filtered-organism%3A%22Homo+sapiens+%28Human%29+%5B9606%5D%22.gff"
# path_uniprot_tab <- "../ref/uniprot/uniprot-filtered-organism%3A%22Homo+sapiens+%28Human%29+%5B9606%5D%22.tab"
# 
# # paths for cached UniprotR objects
# path_taxa_rds   <- file.path(output_dir, "ontology_taxa", "TaxaObj_panelapp.Rds")
# path_go_rds     <- file.path(output_dir, "ontology_taxa", "GeneOntologyObj_panelapp.Rds")
# path_func_rds   <- file.path(output_dir, "ontology_taxa", "ProteinFunction_panelapp.Rds")
# 
# # path for PanelApp genes
# path_panels_rds <- file.path(output_dir, "PanelAppData_genes_combined_Rds")
# 
# # output discussion table
# path_discussion_tsv <- file.path(output_dir, "uniprot_panelapp", "df_uniprot_discussion_panelapp.tsv")
# path_discussion_rds <- file.path(output_dir, "uniprot_panelapp", "df_uniprot_discussion_panelapp.Rds")
# 
# # 1. gene list from PanelAppRex
# 
# df_panels <- readRDS(path_panels_rds)
# 
# gene_list <- df_panels$entity_name |>
#   unique() |>
#   sort()
# 
# cat("PanelApp gene list size:", length(gene_list), "\n")
# 
# # 2. load UniProt meta and tidy gene names
# 
# df_uniprot <- readGFF(path_uniprot_gff)
# 
# df_uniprot_meta <- read.csv(
#   path_uniprot_tab,
#   sep = "\t",
#   check.names = FALSE
# )
# 
# 
# colnames(df_uniprot_meta)[colnames(df_uniprot_meta) == "Entry"] <- "seqid"
# 
# df_uniprot_meta_tidy <- df_uniprot_meta |>
#   tidyr::separate_rows(`Gene names`, sep = " ") |>
#   dplyr::rename(SYMBOL = `Gene names`)
# 
# df_uniprot_meta_tidy <- df_uniprot_meta_tidy |>
#   dplyr::filter(Status == "reviewed") |>
#   dplyr::filter(SYMBOL %in% gene_list)
# 
# cat("Reviewed UniProt rows for PanelApp genes:", nrow(df_uniprot_meta_tidy), "\n")
# 
# Accessions <- df_uniprot_meta_tidy$seqid |> unique()
# 
# cat("Unique UniProt accessions for PanelApp genes:", length(Accessions), "\n")
# 
# # 3. get or load UniprotR objects
# 
# if (file.exists(path_taxa_rds) &&
#     file.exists(path_go_rds) &&
#     file.exists(path_func_rds)) {
#   cat("Loading cached UniprotR objects for PanelApp genes\n")
#   TaxaObj         <- readRDS(path_taxa_rds)
#   GeneOntologyObj <- readRDS(path_go_rds)
#   ProteinFunction <- readRDS(path_func_rds)
# } else {
#   cat("Running UniprotR queries for PanelApp genes\n")
#   TaxaObj         <- GetNamesTaxa(Accessions)
#   GeneOntologyObj <- GetProteinGOInfo(Accessions)
#   ProteinFunction <- GetProteinFunction(Accessions)
#   
#   saveRDS(TaxaObj,         file = path_taxa_rds)
#   saveRDS(GeneOntologyObj, file = path_go_rds)
#   saveRDS(ProteinFunction, file = path_func_rds)
# }
# 
# # 4. merge to build discussion table
# 
# GeneOntologyObj$seqid  <- rownames(GeneOntologyObj)
# TaxaObj$seqid          <- rownames(TaxaObj)
# ProteinFunction$seqid  <- rownames(ProteinFunction)
# 
# data_discussion <- GeneOntologyObj |>
#   dplyr::inner_join(TaxaObj,         by = "seqid") |>
#   dplyr::inner_join(ProteinFunction, by = "seqid")
# 
# df_seqid_gene <- df_uniprot_meta_tidy |>
#   dplyr::select(SYMBOL, seqid) |>
#   dplyr::distinct()
# 
# data_discussion <- df_seqid_gene |>
#   dplyr::inner_join(data_discussion, by = "seqid")
# 
# df_uniprot_discussion <- data_discussion |>
#   dplyr::select(
#     SYMBOL,
#     Protein.names,
#     Gene.Ontology..molecular.function.,
#     Function..CC.
#   ) |>
#   dplyr::distinct()
# 
# cat("Final UniProt discussion table rows:", nrow(df_uniprot_discussion), "\n")
# 
# write_tsv(df_uniprot_discussion, path_discussion_tsv)
# saveRDS(df_uniprot_discussion, path_discussion_rds)
# 
# cat("Saved UniProt discussion data to:\n ")
# cat("TSV:", path_discussion_tsv, "\n ")
# cat("RDS:", path_discussion_rds, "\n")



# build_uniprot_discussion_panelapprex_from_tsv.R

library(dplyr)
library(tidyr)
library(readr)

output_dir <- "../data"
dir.create(file.path(output_dir, "uniprot_panelapp"), showWarnings = FALSE, recursive = TRUE)

path_panels_rds <- file.path(output_dir, "PanelAppData_genes_combined_Rds")

# this is the new TSV you export from UniProt
path_uniprot_tab <- "../ref/uniprot/uniprot_human_reviewed_function.tsv"

path_discussion_tsv <- file.path(output_dir, "uniprot_panelapp", "df_uniprot_discussion_panelapp.tsv")
path_discussion_rds <- file.path(output_dir, "uniprot_panelapp", "df_uniprot_discussion_panelapp.Rds")

# 1. gene list from PanelAppRex
df_panels <- readRDS(path_panels_rds)

gene_list <- df_panels$entity_name |>
  unique() |>
  sort()

cat("PanelApp gene list size:", length(gene_list), "\n")

# 2. load UniProt TSV and tidy gene names

df_uniprot_meta <- read_tsv(
  path_uniprot_tab,
  show_col_types = FALSE
)

names(df_uniprot_meta)
# expect something like:
# "Entry", "Gene Names", "Protein names",
# "Gene Ontology (molecular function)", "Function [CC]", ...

colnames(df_uniprot_meta)[colnames(df_uniprot_meta) == "Entry"] <- "seqid"

df_uniprot_meta_tidy <- df_uniprot_meta |>
  separate_rows(`Gene Names`, sep = " ") |>
  rename(SYMBOL = `Gene Names`)

# 3. restrict to PanelApp genes
df_uniprot_meta_tidy <- df_uniprot_meta_tidy |>
  filter(SYMBOL %in% gene_list)

cat("UniProt rows for PanelApp genes:", nrow(df_uniprot_meta_tidy), "\n")

# 4. build discussion table

# be a bit defensive about column names with spaces
col_prot <- grep("^Protein", names(df_uniprot_meta_tidy), value = TRUE)
col_go_mf <- grep("Gene Ontology \\(molecular function\\)", names(df_uniprot_meta_tidy), value = TRUE)
col_func  <- grep("^Function \\[CC\\]", names(df_uniprot_meta_tidy), value = TRUE)

df_uniprot_discussion <- df_uniprot_meta_tidy |>
  transmute(
    SYMBOL,
    Protein.names = .data[[col_prot]],
    Gene.Ontology..molecular.function. = if (length(col_go_mf) == 1) .data[[col_go_mf]] else NA_character_,
    Function..CC. = if (length(col_func) == 1) .data[[col_func]] else NA_character_
  ) |>
  distinct()

cat("Final UniProt discussion table rows:", nrow(df_uniprot_discussion), "\n")

write_tsv(df_uniprot_discussion, path_discussion_tsv)
saveRDS(df_uniprot_discussion, path_discussion_rds)

cat("Saved UniProt discussion data to:\n ")
cat("TSV:", path_discussion_tsv, "\n ")
cat("RDS:", path_discussion_rds, "\n")