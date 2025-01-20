# Load the required library
library(jsonlite)
library(dplyr)

# Query ----
## Loop pages ----

# Check if the file already exists
if (!file.exists(path_PanelAppData_list_Rds)) {
  # Download and process data if the RDS file does not exist
  
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
  
} else {
  cat("Data files already exist. Skipping download loop.\n")
}

# Load from data to confirm we have it ----
df <- readRDS(file = path_PanelAppData_combined_Rds)

gc()
names(df)

# Get the panel IDs
df_panel_id <- df$id |> 
  sort() |> 
  as.numeric() # currently numbers are used for IDs

print("Line 1 of panel metadata")
head(df, 1)

print("Col names of panel metadata")
names(df)

print("Dimensions of panel metadata")
dim(df)

print("Saved metadata files:")
print(path_PAdata_json)
print(path_PanelAppData_list_Rds)
print(path_PanelAppData_combined_Rds)

