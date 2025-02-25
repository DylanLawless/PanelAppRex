library(shiny)
library(DT)
library(dplyr)

# Rds format
path_data <- "../data"
path_PanelAppData_genes_combined_Rds <- paste0(path_data, "/path_PanelAppData_genes_combined_Rds")
df_core <- readRDS(file= path_PanelAppData_genes_combined_Rds)

df_core <- df_core  |> select(panel_id, entity_name, name, everything()) |> head(1000)

df_small <- df_core %>%
  distinct(panel_id, name)

df_big <- df_core

ui <- fluidPage(
  DTOutput("small_table"),
  DTOutput("big_table")
)

server <- function(input, output, session) {
  output$small_table <- renderDT({
    datatable(
      df_small,
      selection = "multiple",
      options = list(pageLength = 10)
    )
  })
  
  selected_ids <- reactive({
    req(input$small_table_rows_selected)
    df_small$panel_id[input$small_table_rows_selected]
  })
  
  filtered_big <- reactive({
    req(selected_ids())
    df_big %>% filter(panel_id %in% selected_ids())
  })
  
  output$big_table <- renderDT({
    datatable(
      filtered_big(),
      extensions = "Buttons",
      options = list(
        dom = "Bfrtip",
        buttons = list(
          "copy",
          list(
            extend = "csv",
            text = "download filtered csv",
            exportOptions = list(
              modifier = list(
                search = "applied",
                order = "applied"
              )
            )
          )
        ),
        pageLength = 10
      )
    )
  })
}

shinyApp(ui, server)

