library(ggplot2); theme_set(theme_bw())
library(dplyr)
library(tidyr)

# Get your API key token after logging.
# You can login after quick registration.
# Then authorize at: https://panelapp.genomicsengland.co.uk/api/docs/
# The following was produced after testing examples on the API webpage. 

# There seems to be 5 pages and ~1600 ids
# 'https://panelapp.genomicsengland.co.uk/api/v1/panels/?page=%d' 

# Read the CSRF token from the file
api_key_path <- "../credentials/creds_api_key.txt"
csrf_token <- readLines(api_key_path, warn = FALSE)
path_images <-"../images"
path_data <- "../data"
path_PAdata_json <- paste0(path_data, "/PanelAppData.json")
path_PanelAppData_list_Rds <- paste0(path_data, "/PanelAppData_list.Rds")
path_PanelAppData_combined_Rds <- paste0(path_data, "/path_PanelAppData_combined_Rds")

path_PAdata_genes_json <- paste0(path_data, "/PanelAppData_genes.json")
path_PanelAppData_genes_list_Rds <- paste0(path_data, "/PanelAppData_genes_list.Rds")
path_PanelAppData_genes_combined_Rds <- paste0(path_data, "/path_PanelAppData_genes_combined_Rds")

path_PanelAppData_genes_combined_meta <- paste0(path_data, "/PanelAppData_combined_meta")
path_PanelAppData_genes_combined_core <- paste0(path_data, "/PanelAppData_combined_core")
path_PanelAppData_genes_combined_minimal <- paste0(path_data, "/PanelAppData_combined_minimal")



# path_PanelAppData_genes_combined_full <- paste0(path_data, "/PanelAppData_combined_full")

# Run downloads ----
# Get meta
source("1_get_panel_metadata.R")

# select panels to download or skip to get all
# df_panel_id <- head(df_panel_id)

# Get gene panels
source("2_gene_gene_panels.R")

# Prep for analysis ----
df <- readRDS(file = path_PanelAppData_genes_combined_Rds)
df <- df |> select(id, everything())
gc()

# colnames(df)[colnames(df) == 'oldName'] <- 'newName'
colnames(df)[colnames(df) == 'entity_name'] <- 'Gene'
df$id <- as.numeric(df$id)

# save table
df_names <- df %>% select(panel_id, name) %>% unique() # save table. number of panels = name.

df_unique_counts <- df %>%
  summarise_all(~n_distinct(.)) %>%
  pivot_longer(cols = everything(), names_to = "Column_Name", values_to = "Unique_Counts") # save summary of data cols

df_core <- df |>
  select(id,
         mode_of_inheritance,
         Gene,
         confidence_level, 
         name,
         disease_group,
         disease_sub_group,   
         status)

df_minimal <- df |> select(id, Gene)
df_minimal$SYMBOL <- df_minimal$Gene

# export tables ----
write.table(df_names, file= paste0(path_PanelAppData_genes_combined_meta, "_names.tsv"), quote = FALSE, row.names = FALSE, sep = "\t")
write.table(df_unique_counts, file= paste0(path_PanelAppData_genes_combined_meta, "_variable_counts.tsv"), quote = FALSE, row.names = FALSE, sep = "\t")
write.table(df_core, file= paste0(path_PanelAppData_genes_combined_core, ".tsv"), quote = FALSE, row.names = FALSE, sep = "\t")
write.table(df_minimal, file= paste0(path_PanelAppData_genes_combined_minimal, ".tsv"), quote = FALSE, row.names = FALSE, sep = "\t")


# check an example panels
df_select <- df %>% filter(panel_id == 398) # PID genes
df_select <- df %>% filter(panel_id == 467) # Likely inborn error of metabolism
df_select <- df %>% filter(id == 1220) # unexplained child death

# Plots ----
# Count gene per confidence level ----
p1 <- df_core %>%
  select(id, Gene, confidence_level) %>%
  distinct() %>%
  count(confidence_level) %>%
  ggplot(aes(x = confidence_level, y = n) ) +
  geom_bar(stat = "identity", fill = "orange", color = "black") + 
  labs(y = "Number of\ngenes per\nconfidence level")
p1

# No. panels where the same gene is included ----
p2 <- df_core %>%
  select(id, Gene) %>%
  distinct() %>%
  count(Gene) %>%
  arrange(n) %>%
  mutate(gene_index = row_number()) %>%
  ggplot(aes(x = gene_index, y = n)) +
  geom_bar(stat = "identity", fill = "blue") +
  labs(title = "No. panels where the same gene is included",
       x = "Gene index number",
       y = "Number of\npanels with gene")
p2

# No. genes per panel index ----
p3 <- df_core %>%
  select(id, Gene) %>%
  distinct() %>%
  select(id) %>%
  count(id) %>%
  arrange(n) %>%
  mutate(id_index = row_number()) %>%
  ggplot(aes(x = id_index, y = n)) +
  geom_bar(stat = "identity", fill = "purple") +
  labs(title = "Number of genes per panel ID", 
       x = "Panel index number",
       y = "Number of\npanels with gene")

library(patchwork)
patch <- p1 / p2 / p3
patch

ggsave(patch, file = paste0(path_images, "/plot_patch1.pdf") )

# example named
p3_data <- df_core %>%
  select(id, Gene, name) %>%
  distinct() %>%
  select(id, name) %>%
  count(id, name) %>%
  arrange(n) %>%
  mutate(id_index = row_number())
  
p3b <- p3_data %>%
  ggplot(aes(x = id_index, y = n)) +
  geom_bar(stat = "identity", fill = "purple") +
  labs(title = "Number of genes per panel ID", 
       x = "Panel index number",
       y = "Number of panels with gene") +
  ggrepel::geom_label_repel(
    data = filter(p3_data, id == 467 | id == 398),  # Filter df_core to include only ID 467
    aes(label = stringr::str_wrap(name, width = 30)), 
    box.padding = 0.35, point.padding = 0.5, segment.color = 'grey50',
    size = 3, color = "black", nudge_x = -200, nudge_y = 2000
  )

p2_data <- df_core %>%
  select(id, Gene) %>%
  distinct() %>%
  count(Gene) %>%
  arrange(n) %>%
  mutate(gene_index = row_number()) 


p2b <- p2_data %>%
  ggplot(aes(x = gene_index, y = n)) +
  geom_bar(stat = "identity", fill = "blue") +
  labs(title = "No. panels where the same gene is included",
       x = "Gene index number",
       y = "Number of\npanels with gene") + 
ggrepel::geom_label_repel(
  data = filter(p2_data, Gene == "RAG1"),  # Filter df_core to include only ID 467
  aes(label = stringr::str_wrap(Gene, width = 30)), 
  box.padding = 0.35, point.padding = 0.5, segment.color = 'grey50',
  size = 3, color = "black", nudge_x = 0, nudge_y = 20
)

patch2 <- p1 / p2b / p3b
patch2

ggsave(patch2, file = paste0(path_images, "/plot_patch2_annotated_example.pdf") )
ggsave(patch2, file = paste0(path_images, "/plot_patch2_annotated_example.png"), width = 8, height = 6)
# example info ----
print("An example is panel 398 with 572 PID genes which are well established as consensus in the community.")
print("An example is panel 1220 with 1675 genes which are associated with unexplained death in infancy and sudden unexplained death in childhood.")