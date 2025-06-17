library(dplyr)
library(ggplot2)
library(forcats)
library(RColorBrewer)
library(patchwork)
library(stringr)

theme_set(theme_bw())
path_images <- "../images"

# Load data
df_core <- readRDS(file = "../data/path_PanelAppData_genes_combined_Rds")
df_core$Gene <- df_core$entity_name

# Remove NA MOI
# df_core <- df_core %>% filter(!is.na(mode_of_inheritance) & !is.na(Gene))

# Count unique genes
n_genes <- df_core %>% distinct(Gene) %>% nrow()

# Count unique gene + mode_of_inheritance combinations
n_gene_inheritance <- df_core %>% distinct(Gene, mode_of_inheritance) %>% nrow()

# Prepare summary table for mode_of_inheritance
df_plot <- df_core %>%
  distinct(Gene, mode_of_inheritance) %>%
  count(mode_of_inheritance, name = "n") %>%
  mutate(mode_of_inheritance = fct_reorder(mode_of_inheritance, n))

# Plot p1: mode_of_inheritance barplot
p1 <- ggplot(df_plot, aes(x = mode_of_inheritance, y = n, fill = mode_of_inheritance)) +
  geom_col(colour = "black") +
  geom_text(aes(label = n), hjust = -0.2) +
  coord_flip() +
  # scale_fill_brewer(palette = "YlGnBu") +
  scale_fill_brewer(palette = "Spectral") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
  labs(
    x = "mode of inheritance",
    y = "number of unique genes",
    fill = NULL,
    title = NULL
  ) +
  theme(
    legend.position = "none",
    # axis.title.x = element_text(size = 11),
    # axis.title.y = element_text(size = 11)
  )
p1

# Plot p2: total summary bars
df_summary <- tibble::tibble(
  category = c("unique\ngenes", "unique\ngene + MOI"),
  count = c(n_genes, n_gene_inheritance)
)

p2 <- ggplot(df_summary, aes(x = fct_reorder(category, count), y = count, fill = category)) +
  geom_col(colour = "black") +
  geom_text(aes(label = count), vjust = -0.5) +
  scale_fill_brewer(palette = "Dark2") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(
    x = NULL,
    y = "count",
    title = NULL
  ) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 11),
    axis.text.y = element_text(size = 11),
    axis.title.y = element_text(size = 11),
    axis.title.x = element_text(size = 11)
  )
p2

# Combine with patchwork
final_plot <- p2 + p1 + plot_layout(widths = c(1, 1.5)) +
  plot_annotation(
    # title = paste0(
      # "Summary of gene-level inheritance in PanelAppRex base dataset\n",
      # "Unique genes: ", n_genes,
      # " | Gene + inheritance combinations: ", n_gene_inheritance
    # )#,
    # theme = theme(plot.title = element_text(size = 13, face = "bold"))
  ) + plot_annotation(tag_levels = 'A')

print(final_plot)

ggsave(final_plot, file = paste0(path_images, "/plot_patch_uniq_gene_moi_summary.pdf"), width = 8, height = 4)
