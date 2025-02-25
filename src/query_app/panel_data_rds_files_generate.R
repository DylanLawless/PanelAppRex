library(shiny)
library(DT)
library(dplyr)

path_data <- "."
output_directory <- "panel_data_rds_files"
dir.create(file.path(path_data, output_directory), showWarnings = FALSE)

path_PanelAppData_genes_combined_Rds <- paste0(path_data, "/path_PanelAppData_genes_combined_Rds")
df_core <- readRDS(file = path_PanelAppData_genes_combined_Rds)

# Ensure the required columns are selected and retained
df_core <- df_core  |> select(panel_id, entity_name, name, everything())

# Process to save each panel's data into an individual RDS file
unique_panel_ids <- unique(df_core$panel_id)

unique_panel_ids <- unique_panel_ids |> head(20)

for (panel_id in unique_panel_ids) {
  df_subset <- df_core %>% filter(panel_id == panel_id)
  file_name <- paste0(path_data, "/", output_directory, "/path_PanelAppData_genes_combined_Rds_panel_", panel_id, ".Rds")
  saveRDS(df_subset, file = file_name)
}
