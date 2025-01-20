# Part 2 - get panels ----
# Loop ----
# Initialize a list to store data frames for each panel
panel_data_list <- list()

# Loop through panel IDs starting from 1 up to the last_panel_id or until max_panels is reached
# for (panel_id in 1:min(last_panel_id, max_panels)) {
for (panel_id in df_panel_id) {
  
  print(paste("Now running panel ID:", panel_id))
  
  # Construct the curl command with the current panel_id
  curl_command <- sprintf("curl -X 'GET' \\
    'https://panelapp.genomicsengland.co.uk/api/v1/panels/%d/' \\
                        -H 'accept: application/json' \\
                        -H 'X-CSRFToken: %s' > %s", panel_id, csrf_token, path_PAdata_genes_json)
  
  # Execute the curl command to download the JSON file
  system(curl_command)
  
  cat(sprintf("Data for Panel %d saved to temporary slice of %s.\n", panel_id, path_PAdata_genes_json))
  
  # Load the JSON data directly from the file
  panel_data <- fromJSON(path_PAdata_genes_json, flatten = TRUE)
  
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
saveRDS(panel_data_list, file = path_PanelAppData_genes_list_Rds, compress = TRUE)
saveRDS(combined_panel_data, file = path_PanelAppData_genes_combined_Rds, compress = TRUE)

print("Saved panel data files:")
print(path_PAdata_genes_json)
print(path_PanelAppData_genes_list_Rds)
print(path_PanelAppData_genes_combined_Rds)


