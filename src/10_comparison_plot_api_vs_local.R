#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(tidyr)
  library(ggplot2)
  library(scales)
  library(patchwork)
})

# ------------------------------------------------
# Config
# ------------------------------------------------

path_data <- "../data"
path_images <- "../latex/images"
path_plot_pdf <- file.path(path_images, "benchmark_cold_start_manuscript_figure_all_runs.pdf")
path_summary_tsv <- file.path(path_data, "benchmark_cold_start_plot_summary.tsv")

# ------------------------------------------------
# Helpers
# ------------------------------------------------

extract_results_file_meta <- function(path) {
  nm <- basename(path)
  
  m <- str_match(
    nm,
    "^benchmark_results_cold_start_seed([0-9-]+)_plus([0-9]+)_total([0-9]+)_case([0-9]+)\\.tsv$"
  )
  
  if (any(is.na(m))) {
    stop("Could not parse results file name: ", nm, call. = FALSE)
  }
  
  tibble(
    file_path = path,
    seed_panels = m[, 2],
    n_random_extra_panels = as.integer(m[, 3]),
    total_panels = as.integer(m[, 4]),
    case_index = as.integer(m[, 5])
  )
}

clean_method_label <- function(x) {
  x <- as.character(x)
  x <- recode(
    x,
    PanelAppRex_local = "PanelAppRex local",
    GE_PanelApp_API = "GE PanelApp API",
    .default = x
  )
  x
}

method_levels <- c(
  "PanelAppRex local",
  "GE PanelApp API"
)

method_palette <- c(
  "PanelAppRex local" = "#5e3c99",
  "GE PanelApp API" = "#1b9e77"
)

linetype_values <- c(
  "Any returned panel contains causal gene" = "solid",
  "Top-ranked panel contains causal gene" = "22"
)

read_results_files <- function() {
  paths <- list.files(
    path_data,
    pattern = "^benchmark_results_cold_start_seed[0-9-]+_plus[0-9]+_total[0-9]+_case[0-9]+\\.tsv$",
    full.names = TRUE
  )
  
  if (length(paths) == 0) {
    stop("No benchmark_results files found under ", path_data, call. = FALSE)
  }
  
  meta_df <- bind_rows(lapply(paths, extract_results_file_meta))
  
  out <- vector("list", nrow(meta_df))
  for (i in seq_len(nrow(meta_df))) {
    dat <- read_tsv(meta_df$file_path[i], show_col_types = FALSE)
    out[[i]] <- bind_cols(meta_df[i, ], dat)
  }
  
  bind_rows(out)
}

theme_manuscript <- function(base_size = 11) {
  theme_bw(base_size = base_size) +
    theme(
      legend.position = "top",
      legend.title = element_blank(),
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold"),
      strip.background = element_rect(fill = "grey95", colour = "grey75"),
      strip.text = element_text(face = "bold")
    )
}

# ------------------------------------------------
# Load
# ------------------------------------------------

results_all <- read_results_files() |>
  mutate(
    method_plot = clean_method_label(method),
    method_plot = factor(method_plot, levels = method_levels),
    total_panels = as.integer(total_panels),
    case_study = as.character(case_study)
  ) |>
  arrange(total_panels, method_plot, case_study, file_path)

if (!all(method_levels %in% unique(results_all$method_plot))) {
  missing_methods <- setdiff(method_levels, unique(results_all$method_plot))
  warning("Missing expected methods in results: ", paste(missing_methods, collapse = ", "))
}

# ------------------------------------------------
# Summaries
# ------------------------------------------------

summary_by_method_subset <- results_all |>
  group_by(method_plot, total_panels) |>
  summarise(
    n_runs = n(),
    mean_elapsed_total_seconds = mean(elapsed_total_seconds, na.rm = TRUE),
    median_elapsed_total_seconds = median(elapsed_total_seconds, na.rm = TRUE),
    sd_elapsed_total_seconds = sd(elapsed_total_seconds, na.rm = TRUE),
    mean_n_api_calls = mean(n_api_calls, na.rm = TRUE),
    mean_n_candidate_panels_returned = mean(n_candidate_panels_returned, na.rm = TRUE),
    retrieval_success_rate = mean(causal_gene_found_in_any_returned_panel, na.rm = TRUE),
    best_panel_success_rate = mean(best_panel_contains_causal_gene, na.rm = TRUE),
    .groups = "drop"
  ) |>
  arrange(total_panels, method_plot)

# write_tsv(summary_by_method_subset, path_summary_tsv)

plot_df_success <- summary_by_method_subset |>
  select(
    method_plot,
    total_panels,
    retrieval_success_rate,
    best_panel_success_rate
  ) |>
  pivot_longer(
    cols = c(retrieval_success_rate, best_panel_success_rate),
    names_to = "metric",
    values_to = "value"
  ) |>
  mutate(
    metric = recode(
      metric,
      retrieval_success_rate = "Any returned panel contains causal gene",
      best_panel_success_rate = "Top-ranked panel contains causal gene"
    )
  )

# ------------------------------------------------
# Plots
# ------------------------------------------------

p_time <- ggplot() +
  geom_line(
    data = results_all,
    aes(
      x = total_panels,
      y = elapsed_total_seconds,
      group = interaction(method_plot, case_study),
      colour = method_plot
    ),
    alpha = 0.20,
    linewidth = 0.6
  ) +
  geom_point(
    data = results_all,
    aes(
      x = total_panels,
      y = elapsed_total_seconds,
      colour = method_plot
    ),
    alpha = 0.35,
    size = 1.8
  ) +
  geom_line(
    data = summary_by_method_subset,
    aes(
      x = total_panels,
      y = mean_elapsed_total_seconds,
      colour = method_plot
    ),
    linewidth = 1.0
  ) +
  geom_point(
    data = summary_by_method_subset,
    aes(
      x = total_panels,
      y = mean_elapsed_total_seconds,
      colour = method_plot
    ),
    size = 2.8
  ) +
  geom_text(
    data = summary_by_method_subset,
    aes(
      x = total_panels,
      y = mean_elapsed_total_seconds,
      colour = method_plot,
      label = number(mean_elapsed_total_seconds, accuracy = 0.01)
    ),
    vjust = -0.8,
    size = 3,
    show.legend = FALSE
  ) +
  scale_colour_manual(values = method_palette, drop = FALSE) +
  scale_x_continuous(breaks = sort(unique(summary_by_method_subset$total_panels))) +
  scale_y_continuous(labels = label_number(accuracy = 0.01)) +
  labs(
    title = "Cold-start end-to-end time to final result",
    # subtitle = "Faint lines show individual runs. Bold lines show mean by subset size.",
    x = "Number of panels per request",
    y = "Elapsed total\nseconds"
  ) +
  theme_manuscript() +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.3))
  )

p_api <- ggplot(
  summary_by_method_subset,
  aes(x = total_panels, y = mean_n_api_calls, colour = method_plot)
) +
  geom_line(linewidth = 1.0) +
  geom_point(size = 2.8) +
  geom_text(
    aes(label = number(mean_n_api_calls, accuracy = 1)),
    vjust = -0.8,
    size = 3,
    show.legend = FALSE
  ) +
  scale_colour_manual(values = method_palette, drop = FALSE) +
  scale_x_continuous(breaks = sort(unique(summary_by_method_subset$total_panels))) +
  scale_y_continuous(labels = label_number(accuracy = 1)) +
  labs(
    title = "API calls required",
    x = "Number of panels per request",
    y = "Mean API calls"
  ) +
  theme_manuscript() +
  theme(legend.position = "none")+
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.3))
  )

p_candidates <- ggplot(
  summary_by_method_subset,
  aes(x = total_panels, y = mean_n_candidate_panels_returned, colour = method_plot)
) +
  geom_line(linewidth = 1.0) +
  geom_point(size = 2.8) +
  geom_text(
    aes(label = number(mean_n_candidate_panels_returned, accuracy = 0.1)),
    vjust = -0.8,
    size = 3,
    show.legend = FALSE
  ) +
  scale_colour_manual(values = method_palette, drop = FALSE) +
  scale_x_continuous(breaks = sort(unique(summary_by_method_subset$total_panels))) +
  scale_y_continuous(labels = label_number(accuracy = 0.1)) +
  labs(
    title = "Candidate panels returned",
    x = "Number of panels per request",
    y = "Mean candidate\npanels returned"
  ) +
  theme_manuscript() +
  theme(legend.position = "none")+
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.3))
  )

p_success <- ggplot(
  plot_df_success,
  aes(x = total_panels, y = value, colour = method_plot, linetype = metric)
) +
  geom_line(linewidth = 1.0) +
  geom_point(size = 2.4) +
  facet_wrap(~ metric, ncol = 1) +
  scale_colour_manual(values = method_palette, drop = FALSE) +
  scale_linetype_manual(values = linetype_values) +
  scale_x_continuous(breaks = sort(unique(plot_df_success$total_panels))) +
  scale_y_continuous(
    labels = label_percent(accuracy = 1),
    limits = c(0, 1)
  ) +
  labs(
    title = "Retrieval correctness",
    x = "Number of panels per request",
    y = "Rate"
  ) +
  theme_manuscript() +
  theme(
    legend.position = "none"
  )+
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.3))
  )

# combined_plot <- (p_time | p_api) / (p_candidates | p_success) +
combined_plot <- (p_time / p_api / p_candidates) 
  # plot_annotation(
  #   title = "Cold-start benchmark across matched subset sizes",
  #   subtitle = "Primary outcome: time to final ranked result.\nSupporting outcomes: API burden, candidate set size, and retrieval correctness."
  # )

p_time 
p_api
p_candidates
p_success
  
# ------------------------------------------------
# Save
# ------------------------------------------------

ggsave(
  filename = path_plot_pdf,
  plot = combined_plot,
  width = 6,
  height = 8
)

cat("Loaded files:\n")
cat("  results: ", n_distinct(results_all$file_path), "\n", sep = "")
cat("Cases found: ", n_distinct(results_all$case_study), "\n", sep = "")
cat("Subset sizes found: ", paste(sort(unique(results_all$total_panels)), collapse = ", "), "\n", sep = "")
cat("Methods found: ", paste(unique(as.character(results_all$method_plot)), collapse = ", "), "\n", sep = "")
cat("Saved:\n")
cat("  summary: ", path_summary_tsv, "\n", sep = "")
cat("  png: ", path_plot_pdf, "\n", sep = "")
cat("  pdf: ", path_plot_pdf, "\n", sep = "")