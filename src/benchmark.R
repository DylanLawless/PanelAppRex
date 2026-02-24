suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2); theme_set(theme_bw())
  library(patchwork)
  library(RColorBrewer)
  library(stringr)
  library(ggrepel)
})

df <- read.table(
  file = "../data/benchmark.tsv",
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)

total_exact_accuracy  <- mean(df$Causal_gene_in_all_results, na.rm = TRUE) * 100
total_causal_accuracy <- mean(df$Causal_gene_in_subjective_choice, na.rm = TRUE) * 100

cat("Retrieval accuracy for exact_relevance_ratio:", round(total_exact_accuracy, 1), "%\n")
cat("Selection accuracy for Causal_gene_in_subjective_choice:", round(total_causal_accuracy, 1), "%\n")

df_counts <- df[, c(
  "Case_study",
  "Disease_focus",
  "query_result_panels",
  "subjective_best_panels",
  "panels_with_causal_gene"
)]

df_counts$Disease_focus <- factor(
  df_counts$Disease_focus,
  levels = c("Immunology", "Neurology", "Cross-disciplinary")
)

df_long <- pivot_longer(
  df_counts,
  cols = c("query_result_panels",
           "subjective_best_panels",
           "panels_with_causal_gene"),
  names_to = "metric",
  values_to = "value"
)

df_long$metric <- factor(
  df_long$metric,
  levels = c("query_result_panels",
             "panels_with_causal_gene",
             "subjective_best_panels")
)

# annotation only in facet 3
ann_df <- df_long |>
  filter(Disease_focus == "Cross-disciplinary") |>
  summarise(
    Disease_focus = first(Disease_focus),
    Case_study = max(unique(Case_study)),
    .groups = "drop"
  ) |>
  mutate(
    value = 1,
    label = "Successful\ncausal match\nthreshold"
  )

p1 <- ggplot(df_long, aes(x = factor(Case_study), y = value, fill = metric)) +
  geom_bar(
    stat = "identity",
    position = position_dodge(width = 0.9),
    colour = "black"
  ) +
  geom_text(
    aes(label = value),
    position = position_dodge(width = 0.9),
    vjust = -0.5,
    size = 3,
    color = "darkgreen",
  ) +
  geom_hline(
    yintercept = 1,
    colour     = "#88A87B",
    size       = 0.5
  ) +
  geom_hline(
    yintercept = 1,
    linetype   = "dotted",
    colour     = "darkgreen",
    size       = 1
  ) +
  geom_text_repel(
    data = ann_df,
    aes(x = factor(Case_study), y = value, label = label),
    nudge_y = 8,
    segment.colour = "darkgreen",
    colour = "darkgreen",
    size = 3.5,
    inherit.aes = FALSE
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
  labs(
    x = "Case study",
    y = str_wrap("Count", width = 10),
    fill = "Metric"
  ) +
  scale_fill_brewer(
    palette = "RdYlBu",
    labels = function(x) str_wrap(gsub("_", " ", x), width = 40)
  ) +
  facet_wrap(~ Disease_focus, nrow = 1, scales = "free_x") +
  theme_bw() +
  theme(
    legend.position = "top",
    strip.background = element_rect(fill = "white", colour = "black"),
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = 1)
  ) +
  guides(fill = guide_legend(title.position = "top", title = NULL))

p1

# Plots 2 and 3 ----
numeric_cols_percent <- c("Causal_gene_in_all_results", "Causal_gene_in_subjective_choice")

plot_percentages <- function(df, cols) {
  n_cases <- length(unique(df$Case_study))
  pal_15 <- colorRampPalette(RColorBrewer::brewer.pal(9, "YlGnBu"))(max(15, n_cases))
  
  plots <- lapply(cols, function(col) {
    
    annotation <- ""
    if (col == "Causal_gene_in_all_results") {
      annotation <- paste("Retrieval accuracy:", round(total_exact_accuracy, 1), "%")
    } else if (col == "Causal_gene_in_subjective_choice") {
      annotation <- paste("Selection accuracy:", round(total_causal_accuracy, 1), "%")
    }
    
    ggplot(df, aes_string(
      x = "factor(Case_study)",
      y = paste0("(", col, ")*100"),
      fill = "factor(Case_study)"
    )) +
      geom_bar(stat = "identity", colour = "black") +
      geom_text(aes_string(label = paste0("round((", col, ")*100, 0)")), vjust = -0.5, size = 3) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
      labs(
        x = "Case study",
        y = stringr::str_wrap(paste0(gsub("_", " ", col), " (%)"), width = 17)
      ) +
      scale_fill_manual(values = pal_15) +
      guides(fill = "none") +
      annotate("text", x = 1, y = 150, label = annotation, hjust = 0)
  })
  
  plots
}

plots_percent <- plot_percentages(df, numeric_cols_percent)
combined_plot <- patchwork::wrap_plots(plots_percent)
print(combined_plot)

final_plot <- p1 / combined_plot + plot_annotation(tag_levels = 'A')
print(final_plot)

ggsave(final_plot, file = "../latex/images/benchmark_extended.pdf", width = 10, height = 5)
