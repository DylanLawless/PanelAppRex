


library(ggplot2); theme_set(theme_bw())
library(patchwork)
library(tidyr)
library(RColorBrewer)

df <- read.table(file = "../data/benchmark.tsv", header = TRUE, sep = "\t")

total_exact_accuracy <- mean(df$causal_gene_in_all_results) * 100
total_causal_accuracy <- mean(df$causal_gene_in_subjective_choice) * 100

cat("Total accuracy for exact_relevance_ratio:", round(total_exact_accuracy, 1), "%\n")
cat("Total accuracy for causal_gene_in_subjective_choice:", round(total_causal_accuracy, 1), "%\n")

df_counts <- df[, c("Case_study", "query_result_panels", "subjective_best_panels", "panels_with_causal_gene")]

df_long <- pivot_longer(
  df_counts,
  cols = c("query_result_panels", "subjective_best_panels", "panels_with_causal_gene"),
  names_to = "metric",
  values_to = "value"
)

df_long$metric <- factor(df_long$metric, levels = c("query_result_panels", "panels_with_causal_gene", "subjective_best_panels"))

p1 <- ggplot(df_long, aes(x = factor(Case_study), y = value, fill = metric)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9), color = "black") +
  geom_text(aes(label = value),
            position = position_dodge(width = 0.9),
            vjust = -0.5) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
  labs(x = "Case Study", y = stringr::str_wrap("Count", width = 10), fill = "Metric") +
  scale_fill_brewer(palette = "RdYlBu", labels = function(x) stringr::str_wrap(gsub("_", " ", x), width = 40)) +
  theme( legend.position="top") +
  guides(fill=guide_legend(title.position="top"))
    
p1

numeric_cols_percent <- c("causal_gene_in_all_results", "causal_gene_in_subjective_choice")

plot_percentages <- function(df, cols) {
  plots <- lapply(cols, function(col) {
    
    annotation <- ""
    if (col == "causal_gene_in_all_results") {
      annotation <- paste("Total accuracy:", round(total_exact_accuracy, 1), "%")
    } else if (col == "causal_gene_in_subjective_choice") {
      annotation <- paste("Total accuracy:", round(total_causal_accuracy, 1), "%")
    }
    
    ggplot(df, aes_string(x = "factor(Case_study)", y = paste0("(", col, ")*100"), fill = "factor(Case_study)")) +
      geom_bar(stat = "identity", color = "black") +
      geom_text(aes_string(label = paste0("round((", col, ")*100, 0)")), vjust = -0.5) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
      labs(
        x = "Case Study", 
        y = stringr::str_wrap(paste0(gsub("_", " ", col), " (%)"), width = 20)
      ) +
      scale_fill_brewer(palette = "YlGnBu") +
      guides(fill = "none") +
      annotate("text", x = 1, y = 150, label = annotation, hjust = 0)
  })
  return(plots)
}

plots_percent <- plot_percentages(df, numeric_cols_percent)
combined_plot <- wrap_plots(plots_percent)
print(combined_plot)

final_plot <- p1 / combined_plot + plot_annotation(tag_levels = 'A')
print(final_plot)

ggsave(final_plot, file = "../images/benchmark.pdf", width = 8, height = 5)
