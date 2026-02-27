library(dplyr)
library(DT)
library(htmlwidgets)

panelapp_rds_root <- "../data"
output_directory <- "panel_data_html_beta"
dir.create(file.path(panelapp_rds_root, output_directory), showWarnings = FALSE)

# out_dir_rag <- "~/mnt/atlas_data_big/data/panelapprex/uniprot_kb_up000005640_rag"
panelapprex_root  <- "~/mnt/atlas_data_big/data/panelapprex"
out_dir_rag <- file.path(panelapprex_root, "uniprot_kb_up000005640_processed_rag")
# Create a common library folder for all panels
common_lib_dir <- "common_lib"
dir.create(file.path(panelapp_rds_root, output_directory, common_lib_dir), showWarnings = FALSE)

path_PanelAppData_genes_combined_Rds <- file.path(panelapp_rds_root, "PanelAppData_genes_combined_Rds")
df_core <- readRDS(file = path_PanelAppData_genes_combined_Rds)
df_core <- df_core |> select(panel_id, entity_name, name, everything())
df_core$gene_data.gene_name <- gsub(",", " ", df_core$gene_data.gene_name, fixed = TRUE)

# Fix entries with comma separated items
df_core$phenotypes <- vapply(
  df_core$phenotypes,
  function(x) {
    if (length(x) == 0 || all(is.na(x))) return("")
    s <- paste(as.character(x), collapse = ", ")
    s <- gsub("\\s*,\\s*", ", ", s)
    trimws(s)
  },
  character(1)
)

df_core$gene_data.alias <- vapply(
  df_core$gene_data.alias,
  function(x) {
    if (length(x) == 0 || all(is.na(x))) return("")
    s <- paste(as.character(x), collapse = ", ")
    s <- gsub("\\s*,\\s*", ", ", s)
    trimws(s)
  },
  character(1)
)

df_core$publications <- vapply(
  df_core$publications,
  function(x) {
    if (length(x) == 0 || all(is.na(x))) return("")
    x <- trimws(as.character(x))
    x <- x[nzchar(x)]
    if (length(x) == 0) return("")
    paste(unique(x), collapse = " ")
  },
  character(1)
)

# For fast test: keep only panel 398 
# df_core <- df_core |> dplyr::filter(panel_id == 398)
# df_core <- df_core |> dplyr::filter(panel_id == 467)
# df_core <- df_core |> dplyr::filter(panel_id < 20 )

# select cols ----
# cat(names(df_core), sep = "\n")
# df_core <- df_core |> 
#   dplyr::select(
#   name,
#   entity_name,
#   mode_of_inheritance,
#   phenotypes,
#   gene_data.hgnc_symbol,
#   gene_data.hgnc_id,
#   gene_data.omim_gene,
#   gene_data.ensembl_genes.GRch38.90.location,
#   gene_data.ensembl_genes.GRch38.90.ensembl_id,
#   confidence_level,
#   penetrance,
#   gene_data.alias,
#   gene_data.gene_name,
#   gene_data.alias_name,
#   gene_data.gene_symbol,
#   id,
#   disease_group,
#   disease_sub_group,
#   mode_of_inheritance_raw,
#   tags,
#   gene_data.ensembl_genes.GRch37.82.location,
#   gene_data.ensembl_genes.GRch37.82.ensembl_id,
#   publications,
#   panel_id
# )

df_core <- df_core |> 
  dplyr::select(
    entity_name,
    name,
    mode_of_inheritance,
    phenotypes,
    gene_data.hgnc_symbol,
    gene_data.hgnc_id,
    gene_data.omim_gene,
    gene_data.ensembl_genes.GRch38.90.location,
    gene_data.ensembl_genes.GRch38.90.ensembl_id,
    confidence_level,
    penetrance,
    gene_data.alias,
    gene_data.gene_name,
    gene_data.alias_name,
    gene_data.gene_symbol,
    id,
    disease_group,
    disease_sub_group,
    mode_of_inheritance_raw,
    tags,
    gene_data.ensembl_genes.GRch37.82.location,
    gene_data.ensembl_genes.GRch37.82.ensembl_id,
    publications,
    panel_id
  )


# convert to links ----
df_core$gene_data.omim_gene <- vapply(
  df_core$gene_data.omim_gene,
  function(x) {
    
    # Handle NULL / all NA / empty
    if (length(x) == 0 || all(is.na(x))) return("")
    
    x <- paste(x, collapse = ";")
    x <- trimws(x)
    if (!nzchar(x)) return("")
    
    vals <- unique(trimws(unlist(strsplit(x, "[;,]"))))
    vals <- vals[nzchar(vals)]
    if (length(vals) == 0) return("")
    
    links <- vapply(vals, function(v) {
      v_clean <- gsub("[^0-9]", "", v)
      if (!nzchar(v_clean)) return(v)
      sprintf(
        '<a href="https://www.omim.org/entry/%s" target="_blank">%s</a>',
        v_clean, v_clean
      )
    }, character(1))
    
    paste(links, collapse = " ")
  },
  character(1)
)

df_core$gene_data.ensembl_genes.GRch38.90.ensembl_id <- vapply(
  df_core$gene_data.ensembl_genes.GRch38.90.ensembl_id,
  function(x) {
    
    if (length(x) == 0 || all(is.na(x))) return("")
    
    x <- paste(x, collapse = ";")
    x <- trimws(x)
    if (!nzchar(x)) return("")
    
    vals <- unique(trimws(unlist(strsplit(x, "[;,]"))))
    vals <- vals[nzchar(vals)]
    if (length(vals) == 0) return("")
    
    links <- vapply(vals, function(v) {
      v_clean <- trimws(v)
      if (!nzchar(v_clean)) return("")
      sprintf(
        '<a href="http://www.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g=%s" target="_blank">%s</a>',
        v_clean, v_clean
      )
    }, character(1))
    
    paste(links, collapse = " ")
  },
  character(1)
)

df_core$gene_data.hgnc_id <- vapply(
  df_core$gene_data.hgnc_id,
  function(x) {
    
    if (length(x) == 0 || all(is.na(x))) return("")
    
    x <- paste(x, collapse = ";")
    x <- trimws(x)
    if (!nzchar(x)) return("")
    
    vals <- unique(trimws(unlist(strsplit(x, "[;,]"))))
    vals <- vals[nzchar(vals)]
    if (length(vals) == 0) return("")
    
    links <- vapply(vals, function(v) {
      v_clean <- trimws(v)
      if (!nzchar(v_clean)) return("")
      sprintf(
        '<a href="https://www.genenames.org/data/gene-symbol-report/#!/hgnc_id/%s" target="_blank">%s</a>',
        v_clean, v_clean
      )
    }, character(1))
    
    paste(links, collapse = " ")
  },
  character(1)
)

df_core$gene_data.hgnc_symbol <- vapply(
  df_core$gene_data.hgnc_symbol,
  function(x) {
    
    if (length(x) == 0 || all(is.na(x))) return("")
    
    x <- paste(x, collapse = ";")
    x <- trimws(x)
    if (!nzchar(x)) return("")
    
    vals <- unique(trimws(unlist(strsplit(x, "[;,]"))))
    vals <- vals[nzchar(vals)]
    if (length(vals) == 0) return("")
    
    links <- vapply(vals, function(v) {
      v_clean <- trimws(v)
      if (!nzchar(v_clean)) return("")
      sprintf(
        '<a href="https://gnomad.broadinstitute.org/gene/%s" target="_blank">%s</a>',
        v_clean, v_clean
      )
    }, character(1))
    
    paste(links, collapse = " ")
  },
  character(1)
)

# we link to gnomad while showing hgnc gene symbol to save space
names(df_core)[names(df_core) == "gene_data.hgnc_symbol"] <- "hgnc_symbol_gnomad"

# badges ----
badge_gene <- function(x) {
  if (is.na(x) || !nzchar(x)) return("")
  
  vals <- unique(trimws(unlist(strsplit(x, "[;,]"))))
  vals <- vals[nzchar(vals)]
  
  spans <- vapply(vals, function(v) {
    sprintf('<span class="gene-badge">%s</span>', v)
  }, character(1))
  
  paste(spans, collapse = " ")
}

df_core$entity_name <- vapply(
  df_core$entity_name,
  badge_gene,
  character(1)
)

badge_moi <- function(x) {
  if (is.na(x) || !nzchar(x)) return("")
  
  vals <- unique(trimws(unlist(strsplit(x, "[;,]"))))
  vals <- vals[nzchar(vals)]
  
  map_class <- function(v) {
    v2 <- toupper(v)
    if (v2 %in% c("AR", "AUTOSOMAL RECESSIVE")) return("moi-ar")
    if (v2 %in% c("AD", "AUTOSOMAL DOMINANT")) return("moi-ad")
    if (v2 %in% c("XL", "X-LINKED", "X LINKED")) return("moi-xl")
    if (grepl("MITO", v2)) return("moi-mt")
    "moi-other"
  }
  
  spans <- vapply(vals, function(v) {
    cls <- map_class(v)
    sprintf('<span class="moi-badge %s">%s</span>', cls, v)
  }, character(1))
  
  paste(spans, collapse = " ")
}

df_core$mode_of_inheritance <- vapply(
  df_core$mode_of_inheritance,
  badge_moi,
  character(1)
)

unique_panel_ids <- unique(df_core$panel_id) # |> head(20)

# main tables ----
for (pid in unique_panel_ids) {
  # print("making table")
  df_subset <- df_core |> filter(panel_id == pid)
  
  # Create display copy with pretty column names
  df_subset_disp <- df_subset
  names(df_subset_disp) <- gsub("_", " ", names(df_subset_disp), fixed = TRUE)
  names(df_subset_disp) <- gsub("[_.]", " ", names(df_subset_disp))

  
  # dt <- datatable(
  #   df_subset_disp,
  #   escape = FALSE,
  #   extensions = c("Buttons", "FixedColumns"),
  #   options = list(
  #     dom = "Bfrtip",
  #     buttons = c("copy", "csv", "print"),
  #     pageLength = 20,
  #     
  #     scrollX = TRUE,
  #     fixedColumns = list(leftColumns = 2),
  #     
  #     ordering = FALSE,
  #     
  #     language = list(
  #       search = "",
  #       searchPlaceholder = "Enter natural language query..."
  #     ),
  #     initComplete = JS("
  #     function(settings, json) {
  # 
  #       var filter = $('div.dataTables_filter');
  #       filter.css({'width': '100%'});
  #       filter.find('input').css({
  #         'width': '100%',
  #         'border': '1px solid #ccc',
  #         'padding': '8px',
  #         'border-radius': '4px',
  #         'box-shadow': '0 1px 3px rgba(0,0,0,0.2)',
  #         'font-size': '14px',
  #         'font-family': 'system-ui, -apple-system, BlinkMacSystemFont, \"Segoe UI\", Roboto, Oxygen, Ubuntu, Cantarell, \"Open Sans\", \"Helvetica Neue\", sans-serif'
  #       });
  # 
  #       $('table.dataTable').css({
  #         'font-family': 'system-ui, -apple-system, BlinkMacSystemFont, \"Segoe UI\", Roboto, Oxygen, Ubuntu, Cantarell, \"Open Sans\", \"Helvetica Neue\", sans-serif'
  #       });
  # 
  #       if (!document.getElementById('moi-badge-style')) {
  #         var st = document.createElement('style');
  #         st.id = 'moi-badge-style';
  #         st.type = 'text/css';
  #         st.innerHTML = `
  #           .moi-badge {
  #             display: inline-block;
  #             padding: 4px 10px;
  #             margin: 2px 4px 2px 0;
  #             border-radius: 14px;
  #             font-size: 13px;
  #             font-weight: 600;
  #             letter-spacing: 0.2px;
  #             line-height: 1.2;
  #             color: #ffffff;
  #             font-family: system-ui, -apple-system, BlinkMacSystemFont, 'Segoe UI',
  #                          Roboto, Oxygen, Ubuntu, Cantarell, 'Open Sans',
  #                          'Helvetica Neue', sans-serif;
  #           }
  #           .gene-badge {
  #             display: inline-block;
  #             padding: 4px 10px;
  #             margin: 2px 4px 2px 0;
  #             border-radius: 14px;
  #             font-size: 13px;
  #             font-weight: 600;
  #             letter-spacing: 0.2px;
  #             line-height: 1.2;
  #             background: #e3f2fd;
  #             color: #0d47a1;
  #             border: 1px solid rgba(13,71,161,0.15);
  #             font-family: system-ui, -apple-system, BlinkMacSystemFont, 'Segoe UI',
  #                          Roboto, Oxygen, Ubuntu, Cantarell, 'Open Sans',
  #                          'Helvetica Neue', sans-serif;
  #           }
  #           .moi-ar { background: #2e7d32; }
  #           .moi-ad { background: #1565c0; }
  #           .moi-xl { background: #6a1b9a; }
  #           .moi-mt { background: #ef6c00; }
  #           .moi-other { background: #546e7a; }
  # 
  #           .DTFC_LeftWrapper, .DTFC_LeftBodyWrapper, .DTFC_LeftHeadWrapper {
  #             background: #ffffff;
  #           }
  #           .DTFC_LeftBodyWrapper table.dataTable {
  #             background: #ffffff;
  #           }
  #         `;
  #         document.head.appendChild(st);
  #       }
  # 
  #       api.columns.adjust().draw(false);
  #     }
  #   ")
  #   )
  # )
  
  dt <- datatable(
    df_subset_disp,
    escape = FALSE,
    extensions = c("Buttons", "FixedColumns"),
    options = list(
      dom = "Bfrtip",
      buttons = c("copy", "csv", "print"),
      pageLength = 20,
      
      scrollX = TRUE,
      scrollCollapse = TRUE,
      fixedColumns = list(leftColumns = 2),
      
      ordering = FALSE,
      autoWidth = FALSE,
      
      language = list(
        search = "",
        searchPlaceholder = "Enter natural language query..."
      ),
      
      initComplete = JS("
      function(settings, json) {
        var api = new $.fn.dataTable.Api(settings);
        var container = $(api.table().container());

        // Make the whole DataTables widget share one font without overriding other CSS rules
        container.css({
          'font-family': 'system-ui, -apple-system, BlinkMacSystemFont, \"Segoe UI\", Roboto, Oxygen, Ubuntu, Cantarell, \"Open Sans\", \"Helvetica Neue\", sans-serif'
        });

        // Search box styling (scoped to this table only)
        var filter = container.find('div.dataTables_filter');
        filter.css({'width': '100%'});
        filter.find('input').css({
          'width': '100%',
          'border': '1px solid #ccc',
          'padding': '8px',
          'border-radius': '4px',
          'box-shadow': '0 1px 3px rgba(0,0,0,0.2)',
          'font-size': '14px',
          'font-family': 'inherit'
        });

        // Badge styles (your existing block)
        if (!document.getElementById('moi-badge-style')) {
          var st = document.createElement('style');
          st.id = 'moi-badge-style';
          st.type = 'text/css';
          st.innerHTML = `
            .moi-badge {
              display: inline-block;
              padding: 4px 10px;
              margin: 2px 4px 2px 0;
              border-radius: 14px;
              font-size: 13px;
              font-weight: 600;
              letter-spacing: 0.2px;
              line-height: 1.2;
              color: #ffffff;
              font-family: system-ui, -apple-system, BlinkMacSystemFont, 'Segoe UI',
                           Roboto, Oxygen, Ubuntu, Cantarell, 'Open Sans',
                           'Helvetica Neue', sans-serif;
            }
            .gene-badge {
              display: inline-block;
              padding: 4px 10px;
              margin: 2px 4px 2px 0;
              border-radius: 14px;
              font-size: 13px;
              font-weight: 600;
              letter-spacing: 0.2px;
              line-height: 1.2;
              background: #e3f2fd;
              color: #0d47a1;
              border: 1px solid rgba(13,71,161,0.15);
              font-family: system-ui, -apple-system, BlinkMacSystemFont, 'Segoe UI',
                           Roboto, Oxygen, Ubuntu, Cantarell, 'Open Sans',
                           'Helvetica Neue', sans-serif;
            }
            .moi-ar { background: #2e7d32; }
            .moi-ad { background: #1565c0; }
            .moi-xl { background: #6a1b9a; }
            .moi-mt { background: #ef6c00; }
            .moi-other { background: #546e7a; }

            .DTFC_LeftWrapper, .DTFC_LeftBodyWrapper, .DTFC_LeftHeadWrapper {
              background: #ffffff;
            }
            .DTFC_LeftBodyWrapper table.dataTable {
              background: #ffffff;
            }
          `;
          document.head.appendChild(st);
        }

        // Minimal width normalisation for scroll head/body and FixedColumns clones
        if (!document.getElementById('dt-width-patch')) {
          var st2 = document.createElement('style');
          st2.id = 'dt-width-patch';
          st2.type = 'text/css';
          st2.innerHTML = `
            .dataTables_wrapper table.dataTable,
            .dataTables_wrapper table.dataTable thead th,
            .dataTables_wrapper table.dataTable tbody td,
            .dataTables_scrollHead table.dataTable,
            .dataTables_scrollBody table.dataTable,
            .DTFC_LeftHeadWrapper table.dataTable,
            .DTFC_LeftBodyWrapper table.dataTable {
              box-sizing: border-box;
            }

            .dataTables_scrollHeadInner,
            .DTFC_LeftHeadWrapper,
            .DTFC_LeftBodyWrapper {
              box-sizing: border-box;
            }
          `;
          document.head.appendChild(st2);
        }

        // Recalculate widths after styles are in place
        setTimeout(function() {
          api.columns.adjust().draw(false);
          if (api.fixedColumns) { api.fixedColumns().relayout(); }
        }, 0);

        // Keep aligned on resize
        $(window).off('resize.dt-align-' + api.table().node().id);
        $(window).on('resize.dt-align-' + api.table().node().id, function() {
          api.columns.adjust();
          if (api.fixedColumns) { api.fixedColumns().relayout(); }
        });
      }
    ")
 
    )
  )
  file_name_html <- file.path(panelapp_rds_root, output_directory, paste0("panel_", pid, ".html"))
  saveWidget(dt, file = file_name_html, selfcontained = FALSE, libdir = common_lib_dir)
}
