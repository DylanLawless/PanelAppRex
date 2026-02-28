library(biomaRt)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(patchwork)

# Import ----
path_data <- "../data"
path_PanelAppData_genes_combined_Rds <- paste0(path_data, "/PanelAppData_genes_combined_Rds")
df_core <- readRDS(file = path_PanelAppData_genes_combined_Rds)

# Colours ----
cols_update <- RColorBrewer::brewer.pal(n = 10, name = "RdYlBu")[c(5, 1)]

# Helpers ----
fmt_pct <- function(x) {
  out <- sprintf("%.3f", x)
  out <- sub("\\.?0+$", "", out)
  paste0(out, "%")
}

label_count_pct <- function(count, denom) {
  ifelse(
    is.na(denom) | denom == 0,
    as.character(count),
    paste0(count, "\n", fmt_pct(100 * count / denom))
  )
}
# Define denominators ----
df_entries_with_gene <- df_core %>%
  filter(!is.na(entity_name))

genes_unique <- df_entries_with_gene %>%
  distinct(entity_name) %>%
  pull(entity_name)

n_unique_genes <- length(genes_unique)

# Connect to Ensembl once ----
ensembl <- useEnsembl(
  biomart = "genes",
  dataset = "hsapiens_gene_ensembl"
)

# Panel A: HGNC ID (gene-level) ----
hgnc_id_col <- "gene_data.hgnc_id"
has_hgnc_id_col <- hgnc_id_col %in% colnames(df_core)

hgnc_id_before_by_gene <- df_entries_with_gene %>%
  select(entity_name, any_of(hgnc_id_col)) %>%
  distinct(entity_name, .keep_all = TRUE) %>%
  mutate(
    hgnc_id_before = if (has_hgnc_id_col) .data[[hgnc_id_col]] else NA_character_,
    hgnc_id_before = na_if(hgnc_id_before, "")
  ) %>%
  select(entity_name, hgnc_id_before)

n_hgnc_id_before <- sum(!is.na(hgnc_id_before_by_gene$hgnc_id_before))

hgnc_map <- getBM(
  attributes = c("hgnc_symbol", "hgnc_id"),
  filters = "hgnc_symbol",
  values = genes_unique,
  mart = ensembl
) %>%
  distinct(hgnc_symbol, .keep_all = TRUE)

hgnc_id_after_by_gene <- hgnc_id_before_by_gene %>%
  left_join(hgnc_map, by = c("entity_name" = "hgnc_symbol")) %>%
  mutate(
    hgnc_id = na_if(hgnc_id, ""),
    hgnc_id_after = ifelse(!is.na(hgnc_id_before), hgnc_id_before, hgnc_id)
  ) %>%
  select(entity_name, hgnc_id_after)

n_hgnc_id_after <- sum(!is.na(hgnc_id_after_by_gene$hgnc_id_after))

# Plotting table ----
counts_df_hgnc <- tibble(
  metric = factor(
    c("HGNC ID\ntotal entries", "HGNC ID", "Unique genes", "Unique genes"),
    levels = c("HGNC ID\ntotal entries", "Unique genes")
  ),
  stage = factor(
    c("Before", "After", "Before", "After"),
    levels = c("Before", "After")
  ),
  count = c(
    n_hgnc_id_before,
    n_hgnc_id_after,
    n_unique_genes,
    n_unique_genes
  ),
  denom = c(
    n_unique_genes,
    n_unique_genes,
    n_unique_genes,
    n_unique_genes
  )
) %>%
  mutate(
    pct = 100 * count / denom,
    label = label_count_pct(count, denom)
  )

y_max <- max(counts_df_hgnc$count)

# Panel A: HGNC ID (entry-level) + Unique genes (gene-level) ----
hgnc_id_col <- "gene_data.hgnc_id"
has_hgnc_id_col <- hgnc_id_col %in% colnames(df_core)

# Denominators
df_entries_with_gene <- df_core %>%
  filter(!is.na(entity_name))

n_entries_with_genes <- nrow(df_entries_with_gene)

genes_unique <- df_entries_with_gene %>%
  distinct(entity_name) %>%
  pull(entity_name)

n_unique_genes <- length(genes_unique)

# Before: HGNC ID at entry level
hgnc_id_before_by_entry <- df_entries_with_gene %>%
  mutate(
    hgnc_id_before = if (has_hgnc_id_col) .data[[hgnc_id_col]] else NA_character_,
    hgnc_id_before = na_if(hgnc_id_before, "")
  ) %>%
  select(entity_name, hgnc_id_before)

n_hgnc_id_before <- sum(!is.na(hgnc_id_before_by_entry$hgnc_id_before))

# Recover HGNC IDs from HGNC symbol
hgnc_map <- getBM(
  attributes = c("hgnc_symbol", "hgnc_id"),
  filters = "hgnc_symbol",
  values = genes_unique,
  mart = ensembl
) %>%
  distinct(hgnc_symbol, .keep_all = TRUE)

# After: HGNC ID at entry level after fill
hgnc_id_after_by_entry <- hgnc_id_before_by_entry %>%
  left_join(hgnc_map, by = c("entity_name" = "hgnc_symbol")) %>%
  mutate(
    hgnc_id = na_if(hgnc_id, ""),
    hgnc_id_after = ifelse(!is.na(hgnc_id_before), hgnc_id_before, hgnc_id)
  ) %>%
  select(entity_name, hgnc_id_after)

n_hgnc_id_after <- sum(!is.na(hgnc_id_after_by_entry$hgnc_id_after))

# Plotting table
counts_df_hgnc <- tibble(
  metric = factor(
    c("HGNC ID\ntotal entries", "HGNC ID\ntotal entries", "Unique genes", "Unique genes"),
    levels = c("HGNC ID\ntotal entries", "Unique genes")
  ),
  stage = factor(
    c("Before", "After", "Before", "After"),
    levels = c("Before", "After")
  ),
  count = c(
    n_hgnc_id_before,
    n_hgnc_id_after,
    n_unique_genes,
    n_unique_genes
  ),
  denom = c(
    n_entries_with_genes,
    n_entries_with_genes,
    n_unique_genes,
    n_unique_genes
  )
) %>%
  mutate(
    pct = 100 * count / denom,
    label = label_count_pct(count, denom)
  )

y_max <- max(counts_df_hgnc$count)

p1 <- ggplot(counts_df_hgnc, aes(x = metric, y = count, fill = stage)) +
  geom_col(
    position = position_dodge(width = 0.9),
    width = 0.9,
    colour = "black"
  ) +
  geom_text(
    aes(label = label),
    position = position_dodge(width = 0.9),
    vjust = -0.35,
    size = 3.5,
    lineheight = 0.95
  ) +
  scale_y_continuous(
    limits = c(0, y_max * 1.12),
    expand = expansion(mult = c(0, 0.25))
  ) +
  scale_fill_manual(values = cols_update) +
  labs(
    subtitle = "HGNC ID and unique genes",
    x = "",
    y = "\nCount",
    fill = "Stage"
  ) +
  theme_bw() +
  guides(fill = guide_legend(title.position = "top"))

p1

# Panel B: Publications (entry-level) ----

n_entries_with_genes <- df_core %>%
  filter(!is.na(entity_name)) %>%
  nrow()

df_pubs <- df_core %>%
  filter(!is.na(entity_name)) %>%
  mutate(
    has_publication = lengths(publications) > 0
  )

n_entries_with_pubs <- sum(df_pubs$has_publication)

pubs_counts <- tibble(
  metric = factor(
    c("Gene-panel\nentries", "Publication >= 1"),
    levels = c("Gene-panel\nentries", "Publication >= 1")
  ),
  count = c(n_entries_with_genes, n_entries_with_pubs),
  denom = c(n_entries_with_genes, n_entries_with_genes)
) %>%
  mutate(
    pct = 100 * count / denom,
    label = label_count_pct(count, denom)
  )

y_max_pubs <- max(pubs_counts$count)

p2 <- ggplot(pubs_counts, aes(x = metric, y = count, fill = metric)) +
  geom_col(
    width = 0.7,
    colour = "black"
  ) +
  geom_text(
    aes(label = label),
    vjust = -0.35,
    size = 3.5,
    lineheight = 0.95
  ) +
  scale_y_continuous(
    limits = c(0, y_max_pubs * 1.12),
    expand = expansion(mult = c(0, 0.25))
  ) +
  scale_fill_manual(values = c(cols_update[1], cols_update[1])) +
  labs(
    subtitle = "Publication fields",
    x = "",
    y = "\nCount"
  ) +
  theme_bw() +
  theme(legend.position = "none")

p2




# Panel C: Ensembl ID (entry-level) + Unique genes (gene-level) ----
ensembl_id_col <- "gene_data.ensembl_genes.GRch38.90.ensembl_id"
has_ensembl_id_col <- ensembl_id_col %in% colnames(df_core)

# Denominators
df_entries_with_gene <- df_core %>%
  filter(!is.na(entity_name))

n_entries_with_genes <- nrow(df_entries_with_gene)

genes_unique <- df_entries_with_gene %>%
  distinct(entity_name) %>%
  pull(entity_name)

n_unique_genes <- length(genes_unique)

# Before: Ensembl ID at entry level
ensembl_before_by_entry <- df_entries_with_gene %>%
  mutate(
    ensembl_before = if (has_ensembl_id_col) .data[[ensembl_id_col]] else NA_character_,
    ensembl_before = na_if(ensembl_before, "")
  ) %>%
  select(entity_name, ensembl_before)

n_ensembl_before <- sum(!is.na(ensembl_before_by_entry$ensembl_before))

# Recover Ensembl IDs from HGNC symbol
ensembl_map <- getBM(
  attributes = c("hgnc_symbol", "ensembl_gene_id"),
  filters = "hgnc_symbol",
  values = genes_unique,
  mart = ensembl
) %>%
  distinct(hgnc_symbol, .keep_all = TRUE)

# After: Ensembl ID at entry level after fill
ensembl_after_by_entry <- ensembl_before_by_entry %>%
  left_join(ensembl_map, by = c("entity_name" = "hgnc_symbol")) %>%
  mutate(
    ensembl_gene_id = na_if(ensembl_gene_id, ""),
    ensembl_after = ifelse(!is.na(ensembl_before), ensembl_before, ensembl_gene_id)
  ) %>%
  select(entity_name, ensembl_after)

n_ensembl_after <- sum(!is.na(ensembl_after_by_entry$ensembl_after))

# Plotting table
counts_df_ensembl <- tibble(
  metric = factor(
    c("Ensembl ID\ntotal entries", "Ensembl ID\ntotal entries", "Unique genes", "Unique genes"),
    levels = c("Ensembl ID\ntotal entries", "Unique genes")
  ),
  stage = factor(
    c("Before", "After", "Before", "After"),
    levels = c("Before", "After")
  ),
  count = c(
    n_ensembl_before,
    n_ensembl_after,
    n_unique_genes,
    n_unique_genes
  ),
  denom = c(
    n_entries_with_genes,
    n_entries_with_genes,
    n_unique_genes,
    n_unique_genes
  )
) %>%
  mutate(
    pct = 100 * count / denom,
    label = label_count_pct(count, denom)
  )

y_max_ensembl <- max(counts_df_ensembl$count)

p3 <- ggplot(counts_df_ensembl, aes(x = metric, y = count, fill = stage)) +
  geom_col(
    position = position_dodge(width = 0.9),
    width = 0.9,
    colour = "black"
  ) +
  geom_text(
    aes(label = label),
    position = position_dodge(width = 0.9),
    vjust = -0.35,
    size = 3.5,
    lineheight = 0.95
  ) +
  scale_y_continuous(
    limits = c(0, y_max_ensembl * 1.12),
    expand = expansion(mult = c(0, 0.25))
  ) +
  scale_fill_manual(values = cols_update) +
  labs(
    subtitle = "Ensembl ID and unique genes",
    x = "",
    y = "\nCount",
    fill = "Stage"
  ) +
  theme_bw() +
  guides(fill = guide_legend(title.position = "top"))

p3

# Panel D: Disease panel (entry-level) ----

n_entries_with_genes <- df_core %>%
  filter(!is.na(entity_name)) %>%
  nrow()

df_panel_name <- df_core %>%
  filter(!is.na(entity_name)) %>%
  mutate(
    has_panel_name = !is.na(name) & trimws(as.character(name)) != ""
  )

n_entries_with_panel_name <- sum(df_panel_name$has_panel_name)

panel_counts <- tibble(
  metric = factor(
    c("Gene-panel\nentries", "Disease panel"),
    levels = c("Gene-panel\nentries", "Disease panel")
  ),
  count = c(n_entries_with_genes, n_entries_with_panel_name),
  denom = c(n_entries_with_genes, n_entries_with_genes)
) %>%
  mutate(
    pct = 100 * count / denom,
    label = label_count_pct(count, denom)
  )

y_max_panel <- max(panel_counts$count)

p4 <- ggplot(panel_counts, aes(x = metric, y = count, fill = metric)) +
  geom_col(
    width = 0.7,
    colour = "black"
  ) +
  geom_text(
    aes(label = label),
    vjust = -0.35,
    size = 3.5,
    lineheight = 0.95
  ) +
  scale_y_continuous(
    limits = c(0, y_max_panel * 1.12),
    expand = expansion(mult = c(0, 0.25))
  ) +
  scale_fill_manual(values = c(cols_update[1], cols_update[1])) +
  labs(
    subtitle = "Disease panel",
    x = "",
    y = "\nCount"
  ) +
  theme_bw() +
  theme(legend.position = "none")

p4

# Panel E: Mode of inheritance (entry-level) ----

n_entries_with_genes <- df_core %>%
  filter(!is.na(entity_name)) %>%
  nrow()

df_moi <- df_core %>%
  filter(!is.na(entity_name)) %>%
  mutate(
    has_moi = !is.na(mode_of_inheritance) & trimws(as.character(mode_of_inheritance)) != ""
  )

n_entries_with_moi <- sum(df_moi$has_moi)

moi_counts <- tibble(
  metric = factor(
    c("Gene-panel\nentries", "Mode of\ninheritance"),
    levels = c("Gene-panel\nentries", "Mode of\ninheritance")
  ),
  count = c(n_entries_with_genes, n_entries_with_moi),
  denom = c(n_entries_with_genes, n_entries_with_genes)
) %>%
  mutate(
    pct = 100 * count / denom,
    label = label_count_pct(count, denom)
  )

y_max_moi <- max(moi_counts$count)

p5 <- ggplot(moi_counts, aes(x = metric, y = count, fill = metric)) +
  geom_col(
    width = 0.7,
    colour = "black"
  ) +
  geom_text(
    aes(label = label),
    vjust = -0.35,
    size = 3.5,
    lineheight = 0.95
  ) +
  scale_y_continuous(
    limits = c(0, y_max_moi * 1.12),
    expand = expansion(mult = c(0, 0.25))
  ) +
  scale_fill_manual(values = c(cols_update[1], cols_update[1])) +
  labs(
    subtitle = "Mode of inheritance",
    x = "",
    y = "\nCount"
  ) +
  theme_bw() +
  theme(legend.position = "none")

p5

# Panel F: OMIM gene ID (entry-level) ----

n_entries_with_genes <- df_core %>%
  filter(!is.na(entity_name)) %>%
  nrow()

df_omim <- df_core %>%
  filter(!is.na(entity_name)) %>%
  mutate(
    has_omim = !is.na(gene_data.omim_gene) & trimws(as.character(gene_data.omim_gene)) != ""
  )

n_entries_with_omim <- sum(df_omim$has_omim)

omim_counts <- tibble(
  metric = factor(
    c("Gene-panel\nentries", "OMIM gene ID"),
    levels = c("Gene-panel\nentries", "OMIM gene ID")
  ),
  count = c(n_entries_with_genes, n_entries_with_omim),
  denom = c(n_entries_with_genes, n_entries_with_genes)
) %>%
  mutate(
    pct = 100 * count / denom,
    label = label_count_pct(count, denom)
  )

y_max_omim <- max(omim_counts$count)

p6 <- ggplot(omim_counts, aes(x = metric, y = count, fill = metric)) +
  geom_col(
    width = 0.7,
    colour = "black"
  ) +
  geom_text(
    aes(label = label),
    vjust = -0.35,
    size = 3.5,
    lineheight = 0.95
  ) +
  scale_y_continuous(
    limits = c(0, y_max_omim * 1.12),
    expand = expansion(mult = c(0, 0.25))
  ) +
  scale_fill_manual(values = c(cols_update[1], cols_update[1])) +
  labs(
    subtitle = "OMIM gene ID",
    x = "",
    y = "\nCount"
  ) +
  theme_bw() +
  theme(legend.position = "none")

p6

# Joint plot ----
final_plot <- p1 + p3 + p2 + p4 + p5 + p6 +
  # labs(y = "\nCount") +
  plot_layout(guides = "collect", axis = "collect") +
  plot_annotation(tag_levels = "A")  &
  theme(
    axis.text.x = element_text(size = 10)
  )


print(final_plot)

ggsave(final_plot, file = "../latex/images/validation_counts.pdf", width = 12, height = 5)



# 
# New figure legend (Figure S1)
# Validation and recovery of core annotation fields in the PanelAppRex dataset.
# (A) Gene-level completeness of HGNC identifiers, shown as the number of unique genes and the number of unique genes with an HGNC ID, before and after recovery using biomaRt (Ensembl) with HGNC symbol input (entity_name).
# (B) Entry-level completeness of publication annotations, shown as the number of Gene-panel\nentries (rows with a gene) and the number with at least one publication. Not every gene has a publication linked while others have multiple. We do not update with new publications in the core dataset during primary validation.
# (C) Gene-level completeness of Ensembl gene IDs, shown as the number of unique genes and the number of unique genes with an Ensembl gene ID, before and after recovery using biomaRt (Ensembl) with HGNC symbol input. Our primary key is based on gene IDs and a complete match of Ensembl is is not necessary required in the core dataset since we aim for source evidence and We do not update with new data  in the core dataset during primary validation.
# (D) Entry-level completeness of Disease panels.
# (E) Entry-level completeness of mode of inheritance annotations. We do not update with new data  in the core dataset during primary validation.
# (F) Entry-level completeness of OMIM gene identifiers.
# For panels B, D, E, and F, the denominator is the total number of entries with a gene (!is.na(entity_name)). For panels A and C, the denominator is the total number of unique genes. Percentages printed on the relevant bars indicate completeness relative to these denominators.
# Supporting text for Results section
# To verify that PanelAppRex is analysis-ready, we audited completeness of core fields required for downstream joins and filtering (Figure S1). We report entry-level completeness across all gene-panel associations (rows with a gene) for publications, panel names, mode of inheritance, and OMIM gene identifiers, and gene-level completeness for stable identifiers (HGNC ID and Ensembl gene ID) across unique genes.
# Where stable identifiers were missing after source merging, we attempted recovery using Ensembl biomaRt lookups keyed by HGNC symbol. After this recovery step, gene-level completeness of HGNC IDs and Ensembl gene IDs increased to near-complete coverage as shown in Figure S1, while entry-level fields (publications, panel names, mode of inheritance, OMIM) are reported directly from the integrated dataset using explicit denominators.
