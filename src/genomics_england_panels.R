library(ggplot2)

# Get your API key token after logging.
# You can login after quick registration.
# Then authorize at: https://panelapp.genomicsengland.co.uk/api/docs/
# The following was produced after testing examples on the API webpage. 

# There seems to be 5 pages and ~1600 ids
# 'https://panelapp.genomicsengland.co.uk/api/v1/panels/?page=%d' 

# Read the CSRF token from the file
api_key_path <- "../credentials/creds_api_key.txt"
csrf_token <- readLines(api_key_path, warn = FALSE)
path_data <- "../data"
path_PAdata_json <- paste0(path_data, "/PanelAppData.json")
path_PanelAppData_list_Rds <- paste0(path_data, "/PanelAppData_list.Rds")
path_PanelAppData_combined_Rds <- paste0(path_data, "/path_PanelAppData_combined_Rds")

path_PAdata_genes_json <- paste0(path_data, "/PanelAppData_genes.json")
path_PanelAppData_genes_list_Rds <- paste0(path_data, "/PanelAppData_genes_list.Rds")
path_PanelAppData_genes_combined_Rds <- paste0(path_data, "/path_PanelAppData_genes_combined_Rds")


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



# check an example panel
df_select <- df %>%
  filter(panel_id == 398)

# check the count of panels
df_names <- df %>%
  select(panel_id, name) %>% unique()

library(dplyr)
library(tidyr)

# Compute the number of unique values for each column
unique_counts <- df %>%
  summarise_all(~n_distinct(.)) %>%
  pivot_longer(cols = everything(), names_to = "Column_Name", values_to = "Unique_Counts")

# View the dataframe with column names and their unique counts
print(unique_counts)

unique_counts |> filter(Column_Name == "name")


# Plots ----

names(df)

df_sub <- df |>
  select(id,
         entity_name,
         confidence_level, 
         mode_of_inheritance,
         name,
         disease_group,
         disease_sub_group,   
         status)


df_sub %>%
  select(id, entity_name, confidence_level) %>%
  distinct() %>%
  count(confidence_level) %>%
  ggplot(aes(x = confidence_level, y = n) ) +
  geom_bar(stat = "identity") + 
  labs(y = "Number of genes")


df_sub %>%
  select(id, entity_name) %>%
  distinct() %>%
  count(entity_name) %>%
  arrange(n) %>%
  ggplot(aes(x = entity_name, y = n) ) +
  geom_bar(stat = "identity") + 
  labs(title = "Panels per gene", 
       y = "Number of panels with gene")


df_sub$id <- as.character(df_sub$id) 

df_sub %>%
  select(id, entity_name) %>%
  distinct() %>%
  count(id) %>%
  arrange(n) %>%
  ggplot(aes(x = id, y = n) ) +
  geom_bar(stat = "identity") + 
  labs(title = "Genes per panel ID", 
       y = "Number of panels with gene")







# Histogram of Confidence Level Distribution
df_sub %>%
  select(id, confidence_level) %>%
  distinct() %>%
  ggplot(aes(x = confidence_level)) +
  geom_histogram(stat = "count") +
  labs(y = "Number of Genes", x = "Confidence Level")

df_sub %>%
  select(id, entity_name) %>%
  distinct() %>%
  ggplot(aes(x = entity_name)) +
  geom_histogram(stat = "count") +
  labs(y = "Number of Genes", x = "Confidence Level")


df_sub %>%
  # distinct(id) %>%
  ggplot(aes(x = id)) +
  geom_histogram(stat = "count", fill = "blue") +
  theme_minimal()

  

# Histogram of Entity Name Distribution
df_sub %>%
  select(id, entity_name) %>%
  distinct() %>%
  ggplot(aes(x = entity_name)) +
  geom_histogram(stat = "count") +
  labs(title = "Panels per Gene", y = "Number of Panels with Gene", x = "Entity Name")

# Histogram of Gene Occurrence in Panels
df_sub$id <- as.character(df_sub$id)

df_sub %>%
  select(id, entity_name) %>%
  distinct() %>%
  ggplot(aes(x = id)) +
  geom_histogram(stat = "count") +
  labs(title = "Genes per Panel ID", y = "Number of Panels with Gene", x = "Panel ID")








  # scale_fill_manual(values = c("grey", "#ee5d6c"), name = "Carrier\ngenotype", 
                    # guide = guide_legend(reverse = TRUE)) +
  theme_minimal() 
  # xlab("Unique variant\n(arranged by allele count)") +
  # ylab("Allele count") 
  # facet_wrap(~ pathway_id, labeller = labeller(pathway_id = function(x) paste("Pathway ID", x)))


