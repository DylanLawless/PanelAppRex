library(dplyr)
library(DT)
library(htmlwidgets)

# path_root <- "../data"
# path_data <- "."
path_data <- "../data"
output_directory <- "panel_data_html"
dir.create(file.path(path_data, output_directory), showWarnings = FALSE)

# Create a common library folder for all panels
common_lib_dir <- "common_lib"
dir.create(file.path(path_data, output_directory, common_lib_dir), showWarnings = FALSE)

path_PanelAppData_genes_combined_Rds <- file.path(path_data, "PanelAppData_genes_combined_Rds")
df_core <- readRDS(file = path_PanelAppData_genes_combined_Rds)
df_core <- df_core |> select(panel_id, entity_name, name, everything())

unique_panel_ids <- unique(df_core$panel_id)  |> head(20)

for (pid in unique_panel_ids) {
  df_subset <- df_core |> filter(panel_id == pid)
  dt <- datatable(
    df_subset,
    escape = FALSE,
    extensions = "Buttons",
    options = list(
      dom = "Bfrtip",
      buttons = c("copy", "csv", "excel", "pdf", "print"),
      pageLength = 20,
      language = list(
        search = "",
        searchPlaceholder = "Enter natural language query..."
      ),
      initComplete = JS("
        function(settings, json) {
          var filter = $('div.dataTables_filter');
          filter.css({'width': '100%'});
          filter.find('input').css({
            'width': '100%',
            'border': '1px solid #ccc',
            'padding': '8px',
            'border-radius': '4px',
            'box-shadow': '0 1px 3px rgba(0,0,0,0.2)',
            'font-size': '14px',
            'font-family': 'system-ui, -apple-system, BlinkMacSystemFont, \"Segoe UI\", Roboto, Oxygen, Ubuntu, Cantarell, \"Open Sans\", \"Helvetica Neue\", sans-serif'
          });
          $('table.dataTable').css({
            'font-family': 'system-ui, -apple-system, BlinkMacSystemFont, \"Segoe UI\", Roboto, Oxygen, Ubuntu, Cantarell, \"Open Sans\", \"Helvetica Neue\", sans-serif'
          });
        }
      ")
    )
  )
  file_name_html <- file.path(path_data, output_directory, paste0("panel_", pid, ".html"))
  saveWidget(dt, file = file_name_html, selfcontained = FALSE, libdir = common_lib_dir)
}

# the next two code blocks collapse the query data into single cols per panel id row
df_select <- df_core |> 
  select(name, entity_name,
          panel_id, phenotypes, 
         mode_of_inheritance, 
         # mode_of_pathogenicity, 
         disease_group, disease_sub_group,
         gene_data.ensembl_genes.GRch38.90.location,
         gene_data.ensembl_genes.GRch38.90.ensembl_id,
         gene_data.gene_name, gene_data.omim_gene, gene_data.hgnc_id)

df_small <- df_select |>
  group_by(panel_id) |>
  summarise(
    name = first(name),
    gene_count = n_distinct(entity_name),
    entity_names = paste(unique(entity_name), collapse = ";"),
    phenotypes = paste(unique(phenotypes), collapse = ";"),
    mode_of_inheritance = paste(unique(mode_of_inheritance), collapse = ";"),
    # mode_of_pathogenicity = paste(unique(mode_of_pathogenicity), collapse = ";"),
    disease_group = paste(unique(disease_group), collapse = ";"),
    disease_sub_group = paste(unique(disease_sub_group), collapse = ";"),
    gene_data.ensembl_genes.GRch38.90.location = paste(unique(gene_data.ensembl_genes.GRch38.90.location), collapse = ";"),
    gene_data.ensembl_genes.GRch38.90.ensembl_id = paste(unique(gene_data.ensembl_genes.GRch38.90.ensembl_id), collapse = ";"),
    gene_data.gene_name = paste(unique(gene_data.gene_name), collapse = ";"),
    gene_data.omim_gene = paste(unique(gene_data.omim_gene), collapse = ";"),
    gene_data.hgnc_id = paste(unique(gene_data.hgnc_id), collapse = ";"),
  ) |>
  ungroup()

# df_small$link <- sprintf('<a href="./%s/panel_%s.html" target="_blank">open panel</a>',
                         # output_directory, df_small$panel_id)

df_small$name <- sprintf('<a href="./%s/panel_%s.html" target="_blank">%s</a>',
                         output_directory, df_small$panel_id, df_small$name)

df_small <- df_small |> select(name, gene_count, panel_id, everything())

colnames(df_small)[colnames(df_small) == 'gene_count'] <- 'Gene count'
colnames(df_small)[colnames(df_small) == 'name'] <- 'Panel name'

# to html ----
dt_small <- datatable(
  df_small,
  rownames = TRUE,
  escape = FALSE,
  options = list(
    scrollX = F,
    # scrollY = F,  # vertical scroll region
    # scrollY = "200px",  # vertical scroll region
    scroller = TRUE,    # optional for smoother virtual scrolling
    pageLength = 25,
    lengthChange = FALSE,
    deferRender = TRUE,
    autoWidth = FALSE,
    columnDefs = list(
      list(visible = FALSE, targets = seq(3, ncol(df_small)), searchable = TRUE),
      list(width = "20px", targets = 0, orderable = TRUE),
      list(width = "140px", targets = 1),
      list(width = "50px", targets = 2)
    ),
    language = list(
      search = "",
      searchPlaceholder = "Enter natural language query..."
    ),
    initComplete = JS("
      function(settings, json) {
        var filter = $('div.dataTables_filter');
        filter.css({'width': '100%'});
        filter.find('input').css({
          'width': '100%',
          'border': '1px solid #ccc',
          'padding': '8px',
          'border-radius': '4px',
          'box-shadow': '0 1px 3px rgba(0,0,0,0.2)',
          'font-size': '14px',
          'font-family': 'system-ui, -apple-system, BlinkMacSystemFont, \"Segoe UI\", Roboto, Oxygen, Ubuntu, Cantarell, \"Open Sans\", \"Helvetica Neue\", sans-serif'
        });
        $('table.dataTable').css({
          'font-family': 'system-ui, -apple-system, BlinkMacSystemFont, \"Segoe UI\", Roboto, Oxygen, Ubuntu, Cantarell, \"Open Sans\", \"Helvetica Neue\", sans-serif',
          'table-layout': 'fixed'
        });
        $('<style type=\"text/css\"> table.dataTable a { color: #b71c1c; } table.dataTable a:visited { color: #1a237e; } </style>').appendTo('head');
      }
    ")
  )
  
) %>%
  formatStyle(
    'Gene count',
    background = styleColorBar(as.numeric(df_small$`Gene count`), '#43b4eb'),
    backgroundSize = '100% 90%',
    backgroundRepeat = 'no-repeat',
    backgroundPosition = 'center'
  )

dt_small

landing_page <- file.path(paste0(path_data, "/landing_page.html"))
saveWidget(dt_small, landing_page, selfcontained = FALSE)
