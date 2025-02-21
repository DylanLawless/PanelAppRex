# TSV format
path_data <- "../data"
path_PanelAppData_genes_combined_core <- paste0(path_data, "/PanelAppData_combined_core")
path_PanelAppData_genes_combined_minimal <- paste0(path_data, "/PanelAppData_combined_minimal")
df_core <- read.table(file= paste0(path_PanelAppData_genes_combined_core, ".tsv"), sep = "\t")
df_minimal <- read.table(file= paste0(path_PanelAppData_genes_combined_minimal, ".tsv"), sep = "\t")

# Rds format
path_data <- "../data"
path_PanelAppData_genes_combined_Rds <- paste0(path_data, "/path_PanelAppData_genes_combined_Rds")
df_core <- readRDS(file= path_PanelAppData_genes_combined_Rds)
