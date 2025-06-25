library(biomaRt)
library(dplyr)
library(ggplot2)
library(biomaRt)
library(ggplot2)
library(dplyr)

# Import ----
# # TSV format
# path_data <- "../data"
# path_PanelAppData_genes_combined_core <- paste0(path_data, "/PanelAppData_combined_core")
# path_PanelAppData_genes_combined_minimal <- paste0(path_data, "/PanelAppData_combined_minimal")
# df_core <- read.table(file= paste0(path_PanelAppData_genes_combined_core, ".tsv"), sep = "\t")
# df_minimal <- read.table(file= paste0(path_PanelAppData_genes_combined_minimal, ".tsv"), sep = "\t")

# Rds format
path_data <- "../data"
path_PanelAppData_genes_combined_Rds <- paste0(path_data, "/PanelAppData_genes_combined_Rds")
df_core <- readRDS(file= path_PanelAppData_genes_combined_Rds)

# colors ----
# extract two colours from the RdYlBu palette
cols_update <- RColorBrewer::brewer.pal(n = 10, name = "RdYlBu")[c(5, 1)]  # red and blue
cols_unchange <- RColorBrewer::brewer.pal(n = 10, name = "RdYlBu")[c(5, 5)]  # red and blue

# HGNC ----
# counts before
n_genes_before <- n_distinct(df_core$entity_name)
n_hgnc_before <- sum(!is.na(df_core$gene_data.hgnc_symbol))

# identify rows with missing HGNC symbol but present entity_name
df_core_na <- df_core[is.na(df_core$gene_data.hgnc_symbol) & !is.na(df_core$entity_name), ]
genes_to_lookup <- unique(df_core_na$entity_name)

# query Ensembl
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

hgnc_map <- getBM(
  attributes = c("hgnc_symbol", "hgnc_id"),
  filters = "hgnc_symbol",
  values = genes_to_lookup,
  mart = ensembl
)

# merge retrieved HGNCs
df_core_na_updated <- merge(df_core_na, hgnc_map, by.x = "entity_name", by.y = "hgnc_symbol", all.x = TRUE)

# fill in the missing HGNC symbols if recovered
df_core$gene_data.hgnc_symbol[match(df_core_na_updated$entity_name, df_core$entity_name)] <- df_core_na_updated$entity_name

# counts after
n_genes_after <- n_distinct(df_core$entity_name)
n_hgnc_after <- sum(!is.na(df_core$gene_data.hgnc_symbol))

# data for plot
counts_df <- data.frame(
  metric = rep(c("Unique genes", "HGNC ID"), each = 2),
  stage = rep(c("Before", "Updated"), times = 2),
  count = c(n_genes_before, n_genes_after, n_hgnc_before, n_hgnc_after)
)

p1 <- ggplot(counts_df, aes(x = metric, y = count, fill = stage)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9), colour = "black") +
  geom_text(aes(label = count), 
            position = position_dodge(width = 0.9), 
            vjust = -0.5, size = 3.5) +
  scale_y_continuous(limits =  c(0,69000)) +
  scale_fill_manual(values = cols_update) +
  labs(
    subtitle = "Gene and HGNC ID",
    x = "",
    y = "Count",
    fill = "Stage"
  ) +
  theme_bw() +
  guides(fill = guide_legend(title.position = "top"))

p1

# Publications -----

# count entries with non-missing genes
n_entries_with_genes <- df_core %>%
  filter(!is.na(entity_name)) %>%
  nrow()

# count entries with non-missing publications
n_entries_with_pubs <- df_core %>%
  filter(!is.na(publications)) %>%
  nrow()

# prepare data
pubs_counts <- data.frame(
  metric = c("Gene entry", "Publication\n>= 1"),
  count = c(n_entries_with_genes, n_entries_with_pubs)
)

# plot
p2 <- ggplot(pubs_counts, aes(x = metric, y = count, fill = metric)) +
  geom_bar(stat = "identity", colour = "black", width = 0.7) +
  geom_text(aes(label = count), vjust = -0.5, size = 3.5) +
  scale_y_continuous(limits =  c(0,69000)) +
    # scale_fill_manual(values = cols) +
  scale_fill_manual(values = cols_unchange) +
  labs(
    subtitle = "Publication fields",
    x = "",
    y = "Count",
    fill = "Field"
  ) +
  theme_bw() +
  theme(legend.position = "none")

p2

# Ensembl ----
# Count before update
n_genes_before <- n_distinct(df_core$entity_name)
n_ensembl_before <- sum(!is.na(df_core$gene_data.ensembl_genes.GRch38.90.ensembl_id))

# Identify missing Ensembl IDs with known gene names
df_core_missing_ensembl <- df_core %>%
  filter(is.na(gene_data.ensembl_genes.GRch38.90.ensembl_id) & !is.na(entity_name))

genes_to_lookup <- unique(df_core_missing_ensembl$entity_name)

# Query Ensembl via biomaRt
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

ensembl_map <- getBM(
  attributes = c("hgnc_symbol", "ensembl_gene_id"),
  filters = "hgnc_symbol",
  values = genes_to_lookup,
  mart = ensembl
)

# Merge and fill missing Ensembl IDs
df_core_missing_ensembl_updated <- merge(
  df_core_missing_ensembl,
  ensembl_map,
  by.x = "entity_name",
  by.y = "hgnc_symbol",
  all.x = TRUE
)

match_idx <- match(df_core_missing_ensembl_updated$entity_name, df_core$entity_name)
df_core$gene_data.ensembl_genes.GRch38.90.ensembl_id[match_idx] <- df_core_missing_ensembl_updated$ensembl_gene_id

# Count after update
n_ensembl_after <- sum(!is.na(df_core$gene_data.ensembl_genes.GRch38.90.ensembl_id))

# Prepare counts for plotting
counts_df_ensembl <- data.frame(
  metric = rep(c("Unique genes", "Ensembl ID"), each = 2),
  stage = rep(c("Before", "Updated"), times = 2),
  count = c(n_genes_before, n_genes_before, n_ensembl_before, n_ensembl_after)
)

# Plot
p3 <- ggplot(counts_df_ensembl, aes(x = metric, y = count, fill = stage)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9), colour = "black") +
  geom_text(aes(label = count),
            position = position_dodge(width = 0.9),
            vjust = -0.5, size = 3.5) +
  scale_y_continuous(limits =  c(0,69000)) +
  scale_fill_manual(values = cols_update) +
  labs(
    subtitle = "Gene and Ensembl ID",
    x = "",
    y = "Count",
    fill = "Stage"
  ) +
  theme_bw() +
  guides(fill = guide_legend(title.position = "top"))

p3

# Disease panel name ----
n_entries_with_panel_name <- df_core %>%
  filter(!is.na(name)) %>%
  nrow()

panel_counts <- data.frame(
  metric = c("Gene entry", "Panel name"),
  count = c(n_entries_with_genes, n_entries_with_panel_name)
)

p4 <- ggplot(panel_counts, aes(x = metric, y = count, fill = metric)) +
  geom_bar(stat = "identity", colour = "black", width = 0.7) +
  geom_text(aes(label = count), vjust = -0.5, size = 3.5) +
  scale_y_continuous(limits =  c(0,69000)) +
  scale_fill_manual(values = cols_unchange) +
  labs(
    subtitle = "Gene and disease\npanel name",
    x = "",
    y = "Count",
    fill = "Field"
  ) +
  theme_bw() +
  theme(legend.position = "none")

p4

# Disease MOI ----
n_entries_with_moi <- df_core %>%
  filter(!is.na(mode_of_inheritance)) %>%
  nrow()

moi_counts <- data.frame(
  metric = c("Gene entry", "Mode of\ninheritance"),
  count = c(n_entries_with_genes, n_entries_with_moi)
)

p5 <- ggplot(moi_counts, aes(x = metric, y = count, fill = metric)) +
  geom_bar(stat = "identity", colour = "black", width = 0.7) +
  geom_text(aes(label = count), vjust = -0.5, size = 3.5) +
  scale_y_continuous(limits =  c(0,69000)) +
  scale_fill_manual(values = cols_unchange) +
  labs(
    subtitle = "Mode of inheritance",
    x = "",
    y = "Count",
    fill = "Field"
  ) +
  theme_bw() +
  theme(legend.position = "none")

p5

# OMIM ----

# OMIM gene ----
n_entries_with_omim <- df_core %>%
  filter(!is.na(gene_data.omim_gene)) %>%
  nrow()

omim_counts <- data.frame(
  metric = c("Gene entry", "OMIM ID"),
  count = c(n_entries_with_genes, n_entries_with_omim)
)

p6 <- ggplot(omim_counts, aes(x = metric, y = count, fill = metric)) +
  geom_bar(stat = "identity", colour = "black", width = 0.7) +
  geom_text(aes(label = count), vjust = -0.5, size = 3.5) +
  scale_y_continuous(limits = c(0, 69000)) +
  scale_fill_manual(values = cols_unchange) +
  labs(
    subtitle = "Gene and OMIM ID",
    x = "",
    y = "Count",
    fill = "Field"
  ) +
  theme_bw() +
  theme(legend.position = "none")

p6



# annotate ----
add_100pct <- function(p) {
  ymax <- ggplot_build(p)$layout$panel_params[[1]]$y.range[2]
  y_annot <- ymax * 0.95
  p + annotate("text", x = Inf, y = y_annot, label = "100%", hjust = 1.1, vjust = 1.1)
}

# apply to each plot
p1 <- add_100pct(p1)
p2 <- add_100pct(p2)
p3 <- add_100pct(p3)
p4 <- add_100pct(p4)
p5 <- add_100pct(p5)
p6 <- add_100pct(p6)

# joint plot ----
library(patchwork)
final_plot <- p1 + p2 + p3 + p4 + p5 + p6 + plot_layout(guides = 'collect', axis = "collect")  + plot_annotation(tag_levels = 'A')

print(final_plot)

ggsave(final_plot, file = "../images/validation_counts.pdf", width = 11, height = 5.5)



