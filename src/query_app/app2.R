library(shiny)
library(DT)
library(dplyr)

output_directory <- "panel_data_rds_files"
path_data <- "."
path_PanelAppData_genes_combined_Rds <- paste0(path_data, "/path_PanelAppData_genes_combined_Rds")
df_core <- readRDS(file= path_PanelAppData_genes_combined_Rds)
df_core <- df_core  |> select(panel_id, entity_name, name, everything()) # |> head(1000)


df_small <- df_core %>%
  group_by(panel_id) %>%
  summarise(
    name = first(name),
    entity_names = paste(unique(entity_name), collapse = ", "),
    gene_count = n_distinct(entity_name)
  ) %>%
  ungroup()

# Define the user interface
ui <- fluidPage(
    mainPanel(
      DTOutput("small_table"),
      DTOutput("big_table")
    )
  )


# Define server logic
server <- function(input, output, session) {
  # Render the small table
  output$small_table <- renderDT({
    datatable(
      df_small,
      selection = "single",
      options = list(pageLength = 10)
    )
  })
  
  # Reactive value to hold the selected panel ID
  selected_id <- reactive({
    input$small_table_rows_selected
  })
  
  # Load the corresponding panel data when a row is selected
  output$big_table <- renderDT({
    req(selected_id())  # Ensure there is a selection
    panel_id <- df_small$panel_id[selected_id()]
    file_path <- paste0(path_data, "/", output_directory, "/path_PanelAppData_genes_combined_Rds_panel_", panel_id, ".Rds")
    df_panel <- readRDS(file = file_path)
    
    datatable(
      df_panel,
      extensions = "Buttons",
      options = list(
        dom = 'Bfrtip',
        buttons = c('copy', 'csv', 'excel', 'pdf', 'print'),
        pageLength = 10
      )
    )
  })
}

# Run the application
shinyApp(ui, server)
