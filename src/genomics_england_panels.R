# Load the required library
library(jsonlite)
library(dplyr)

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

# Query ----
## Loop pages ----

# Instead of finding panel lets get all pages and cross check

# Initialize a list to store data frames for each page
page_data_list <- list()

# Define a maximum number of iterations to prevent endless loops
max_page <- 5
last_page_id <- 6 # current number of panels reported

# Loop through panel IDs starting from 1 up to the last_page_id or until max_page is reached
for (page_id in 1:min(last_page_id, max_page)) {
  
  # Correctly format the curl command using %s for string substitution
  # curl_command <- sprintf("curl -X 'GET' \\
  # 'https://panelapp.genomicsengland.co.uk/api/v1/panels/?page=%d' \\
  # -H 'accept: application/json' \\
  # -H 'X-CSRFToken: %s' > PanelAppData.json", page_id, csrf_token)
  
  # Correctly format the curl command using %s for string substitution, including the output file path
  curl_command <- sprintf("curl -X 'GET' \\
  'https://panelapp.genomicsengland.co.uk/api/v1/panels/?page=%d' \\
                        -H 'accept: application/json' \\
                        -H 'X-CSRFToken: %s' > %s", page_id, csrf_token, path_PAdata_json)
  
  system(curl_command)
  cat("Data saved to tempory slice of:", path_PAdata_json)
  
  # Load the JSON data directly from the file
  panel_data <- fromJSON(path_PAdata_json, flatten = TRUE)
  panel_data$results$name
  panel_data$results
  # Extract results - no gene data - into a DataFrame
  gene_df <- as.data.frame(panel_data$results)
  panel_data$count
  panel_data$`next`
  panel_data$previous
  
  # Add the page_id as a column to distinguish data from different panels
  gene_df$page_id <- page_id
  
  # Store the dataframe in the list with page_id as the name
  page_data_list[[as.character(page_id)]] <- gene_df
  
  
  # Break the loop if the last intended panel ID is reached
  if (page_id == last_page_id) {
    break
  }
}

# Combine all dataframes into a single dataframe if any data exists in the list
if (length(page_data_list) > 0) {
  combined_page_data <- bind_rows(page_data_list)
  # Check the combined dataframe
  print(head(combined_page_data))
} else {
  print("No data was combined; all panels fetched were empty.")
}

# Store copies ----
saveRDS(page_data_list, file = path_PanelAppData_list_Rds, compress = TRUE)
saveRDS(combined_page_data, file = path_PanelAppData_combined_Rds, compress = TRUE)
rm(combined_page_data, page_data_list)

# Load from data to confirm we have it ----
df <- readRDS(file = path_PanelAppData_combined_Rds)

gc()
names(df)





















# Assuming you want to test the command for a specific panel ID, e.g., 1
panel_id <- 1

# Correctly format the curl command using %s for string substitution
curl_command <- sprintf("curl -X 'GET' \\
  'https://panelapp.genomicsengland.co.uk/api/v1/panels/%d/' \\
                        -H 'accept: application/json' \\
                        -H 'X-CSRFToken: %s' > PanelAppData.json", panel_id, csrf_token)

system(curl_command)
cat("Data saved to 'PanelAppData.json'.")

# Load the JSON data directly from the file
panel_data <- fromJSON("PanelAppData.json", flatten = TRUE)

# Extract the gene data into a DataFrame
gene_df <- as.data.frame(panel_data$genes)

# Define the list of attributes you want to add to the gene DataFrame
attributes_to_add <- c("id", "hash_id", "name", "disease_group", "disease_sub_group", "status", "version", "version_created")

# the following attributes sometimes are lists and are dropped here for simplicity:
# "stats", "types", "genes", 
# "strs", "regions", "relevant_disorders"

# Loop over the list of attributes and add each one as a new column
for (attr in attributes_to_add) {
  # Check if the attribute exists in panel_data
  if (attr %in% names(panel_data)) {
    # Add the attribute as a new column in gene_df
    gene_df[[attr]] <- ifelse(is.null(panel_data[[attr]]), NA, panel_data[[attr]])
  }
}

# View the DataFrame to confirm additions
print(head(gene_df))
rm(gene_df, panel_data)

# Loop ----
# Initialize a list to store data frames for each panel
panel_data_list <- list()

# Define a maximum number of iterations to prevent endless loops
max_panels <- 3000
last_panel_id <- 451 # current number of panels reported

# Loop through panel IDs starting from 1 up to the last_panel_id or until max_panels is reached
for (panel_id in 1:min(last_panel_id, max_panels)) {
  
  # Construct the curl command with the current panel_id
  curl_command <- sprintf("curl -X 'GET' \\
  'https://panelapp.genomicsengland.co.uk/api/v1/panels/%d/' \\
                        -H 'accept: application/json' \\
                        -H 'X-CSRFToken: %s' > temp_PanelAppData.json", panel_id, csrf_token)
  
  # Execute the curl command to download the JSON file
  system(curl_command)
  cat(sprintf("Data for Panel %d saved to 'temp_PanelAppData.json'.\n", panel_id))
  
  # Load the JSON data directly from the file
  panel_data <- fromJSON("temp_PanelAppData.json", flatten = TRUE)
  
  # Check if the genes data is available and not empty
  if (!is.null(panel_data$genes) && length(panel_data$genes) > 0) {
    # Extract the gene data into a DataFrame
    gene_df <- as.data.frame(panel_data$genes)
    
    # Define the list of attributes you want to add to the gene DataFrame
    attributes_to_add <- c("id", "hash_id", "name", "disease_group", "disease_sub_group", "status", "version", "version_created")

        # the following attributes sometimes are lists and are dropped here for simplicity:
    # "stats", "types", "genes", 
    # "strs", "regions", "relevant_disorders"

    # Loop over the list of attributes and add each one as a new column
    for (attr in attributes_to_add) {
      if (attr %in% names(panel_data)) {
        gene_df[[attr]] <- ifelse(is.null(panel_data[[attr]]), NA, panel_data[[attr]])
      }
    }
    
    # Add the panel_id as a column to distinguish data from different panels
    gene_df$panel_id <- panel_id
    
    # Store the dataframe in the list with panel_id as the name
    panel_data_list[[as.character(panel_id)]] <- gene_df
  } else {
    cat(sprintf("No data found for Panel %d or data is incomplete.\n", panel_id))
  }
  
  # Break the loop if the last intended panel ID is reached
  if (panel_id == last_panel_id) {
    break
  }
}

# Combine all dataframes into a single dataframe if any data exists in the list
if (length(panel_data_list) > 0) {
  combined_panel_data <- bind_rows(panel_data_list)
  # Check the combined dataframe
  print(head(combined_panel_data))
} else {
  print("No data was combined; all panels fetched were empty.")
}

# Combine all dataframes into a single dataframe if any data exists in the list
if (length(panel_data_list) > 0) {
  combined_panel_data <- bind_rows(panel_data_list)
  # Check the combined dataframe
  print(head(combined_panel_data))
} else {
  print("No data was combined; all panels fetched were empty.")
}

# Store copies ----
saveRDS(panel_data_list, file = "../data/panel_data_list.Rds", compress = TRUE)
saveRDS(combined_panel_data, file = "../data/panel_data_combined.Rds", compress = TRUE)

# Prep for analysis ----
combined_panel_data <- readRDS(file = "../data/panel_data_combined.Rds")
df <- combined_panel_data
rm(combined_panel_data)
gc()
max(df$panel_id)

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

# # Select specific columns to view
# df <- df %>% 
#   select(entity_name, 
#          gene_symbol = gene_data.gene_symbol, 
#          ensembl_id = gene_data.ensembl_genes.GRch38.90.ensembl_id, everything())
# 
# # Print the selected information
# head(df)


# Loop for pages instead. 

# Assuming you want to test the command for a specific panel ID, e.g., 1
# page_id <- 5



