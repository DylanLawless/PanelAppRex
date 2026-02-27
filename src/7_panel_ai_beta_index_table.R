library(dplyr)
library(DT)
library(htmlwidgets)
library(jsonlite)

panelapp_rds_root <- "../data"
output_directory <- "panel_data_html_beta"
dir.create(file.path(panelapp_rds_root, output_directory), showWarnings = FALSE)

panelapprex_root  <- "~/mnt/atlas_data_big/data/panelapprex"
out_dir_rag <- file.path(panelapprex_root, "uniprot_kb_up000005640_processed_rag")

common_lib_dir <- "common_lib"
dir.create(file.path(panelapp_rds_root, output_directory, common_lib_dir), showWarnings = FALSE)

path_PanelAppData_genes_combined_Rds <- file.path(panelapp_rds_root, "PanelAppData_genes_combined_Rds")
df_core <- readRDS(file = path_PanelAppData_genes_combined_Rds)
df_core <- df_core |> select(panel_id, entity_name, name, everything())

# Apply same optional filters if desired
# df_core <- df_core |> dplyr::filter(panel_id == 398)
# df_core <- df_core |> dplyr::filter(panel_id < 20)

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
    name = dplyr::first(name),
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

df_small$name <- sprintf('<a href="./%s/panel_%s.html" target="_blank">%s</a>',
                         output_directory, df_small$panel_id, df_small$name)

df_small <- df_small |> select(name, gene_count, panel_id, everything())

colnames(df_small)[colnames(df_small) == 'gene_count'] <- 'Gene count'
colnames(df_small)[colnames(df_small) == 'name'] <- 'Panel name'


# make background transparent  ----

# Beta UI add on: compact clickable info icon (popover) for AI RAG text ----
rag_read_txt_keep_newlines <- function(path) {
  if (!file.exists(path)) return(NA_character_)
  txt <- readLines(path, warn = FALSE, encoding = "UTF-8")
  txt <- paste(txt, collapse = "\n")
  if (!nzchar(trimws(txt))) return(NA_character_)
  txt
}

rag_over_path <- function(pid) file.path(out_dir_rag, sprintf("panel_%d_rag_overview.txt", pid))
rag_res_path  <- function(pid) file.path(out_dir_rag, sprintf("panel_%d_rag_result.txt", pid))

pid_vec <- df_small$panel_id

rag_over_raw <- vapply(
  pid_vec,
  function(pid) rag_read_txt_keep_newlines(rag_over_path(pid)),
  character(1)
)
rag_res_raw  <- vapply(
  pid_vec,
  function(pid) rag_read_txt_keep_newlines(rag_res_path(pid)),
  character(1)
)

format_bullets_html <- function(x) {
  if (is.na(x) || !nzchar(trimws(x))) return("")
  lines <- unlist(strsplit(x, "\n", fixed = TRUE))
  lines <- trimws(lines)
  lines <- lines[nzchar(lines)]
  lines <- sub("^[-•\\*]\\s*", "", lines)
  if (length(lines) == 0) return("")
  paste0("<ul class=\"rag-ul\"><li>", paste(lines, collapse = "</li><li>"), "</li></ul>")
}

rag_res_html <- vapply(rag_res_raw, format_bullets_html, character(1))

df_small2 <- df_small
df_small2$rag_overview_raw <- rag_over_raw
df_small2$rag_result_html  <- rag_res_html

df_small2$rag_info <- ifelse(
  !is.na(df_small2$rag_overview_raw) & nzchar(trimws(df_small2$rag_overview_raw)),
  sprintf(
    '<span class="rag-info-icon" data-pid="%s" role="button" tabindex="0" aria-label="Show AI summary">&#9432;</span>',
    df_small2$panel_id
  ),
  ""
)

colnames(df_small2)[colnames(df_small2) == "rag_info"] <- "Info"

df_small2 <- df_small2 |>
  select(
    `Panel name`,
    `Gene count`,
    panel_id,
    Info,
    everything()
  )

rag_payload <- df_small2 |>
  transmute(
    panel_id,
    rag_overview_raw = ifelse(is.na(rag_overview_raw), "", rag_overview_raw),
    rag_result_html  = ifelse(is.na(rag_result_html), "", rag_result_html)
  ) |>
  jsonlite::toJSON(auto_unbox = TRUE)

popover_js <- JS(sprintf("
  (function() {
    var ragData = %s;
    var pop = null;
    var popPid = null;

    function byPid(pid) {
      for (var i = 0; i < ragData.length; i++) {
        if (String(ragData[i].panel_id) === String(pid)) return ragData[i];
      }
      return null;
    }

    function closePop() {
      if (pop) {
        pop.remove();
        pop = null;
        popPid = null;
      }
    }

    function clamp(val, lo, hi) {
      return Math.max(lo, Math.min(hi, val));
    }

    function installStyles() {
      if (document.getElementById('rag-popover-style')) return;
      var st = document.createElement('style');
      st.id = 'rag-popover-style';
      st.type = 'text/css';
      st.textContent = `
  .rag-info-icon {
    cursor: pointer;
    color: #1565c0;
    font-size: 14px;
    line-height: 14px;
    display: inline-block;
    width: 14px;
    text-align: center;
    user-select: none;
    font-family: system-ui, -apple-system, BlinkMacSystemFont, 'Segoe UI',
                 Roboto, Oxygen, Ubuntu, Cantarell, 'Open Sans',
                 'Helvetica Neue', sans-serif;
  }
  .rag-popover {
    position: absolute;
    z-index: 20000;
    background: #ffffff;
    border: 1px solid rgba(0,0,0,0.18);
    border-radius: 6px;
    box-shadow: 0 4px 14px rgba(0,0,0,0.18);
    padding: 10px 12px;
    width: 50vw;
    max-width: calc(100vw - 24px);
    max-height: 60vh;
    overflow: auto;
    font-family: system-ui, -apple-system, BlinkMacSystemFont, 'Segoe UI',
                 Roboto, Oxygen, Ubuntu, Cantarell, 'Open Sans',
                 'Helvetica Neue', sans-serif;
    font-size: 14px;
    line-height: 1.4;
  }
  .rag-h {
    font-weight: 600;
    margin: 2px 0 6px 0;
    font-size: 14px;
  }
  .rag-p {
    margin: 0 0 10px 0;
    font-size: 14px;
  }
  .rag-ul {
    margin: 6px 0 10px 18px;
    padding: 0;
  }
  .rag-ul li {
    margin: 2px 0;
    font-size: 14px;
  }
  .rag-close {
    font-size: 12px;
    color: #c62828;
    cursor: pointer;
    text-align: right;
    margin-top: 6px;
    user-select: none;
    font-family: system-ui, -apple-system, BlinkMacSystemFont, 'Segoe UI',
                 Roboto, Oxygen, Ubuntu, Cantarell, 'Open Sans',
                 'Helvetica Neue', sans-serif;
  }
`;
      document.head.appendChild(st);
    }

    function openPop(anchorEl, pid) {
      if (pop && popPid === pid) {
        closePop();
        return;
      }

      closePop();
      var row = byPid(pid);
      if (!row) return;

      pop = document.createElement('div');
      pop.className = 'rag-popover';
      pop.setAttribute('role', 'dialog');
      pop.setAttribute('aria-label', 'AI summary');

      var html = '';

      if (row.rag_overview_raw && row.rag_overview_raw.trim().length > 0) {
        var safeOverview = row.rag_overview_raw
          .replace(/&/g, '&amp;')
          .replace(/</g, '&lt;')
          .replace(/>/g, '&gt;')
          .replace(/\\n/g, '<br>');
        html += '<div class=\"rag-h\">Overview</div><div class=\"rag-p\">' + safeOverview + '</div>';
      }

      if (row.rag_result_html && row.rag_result_html.trim().length > 0) {
        html += '<div class=\"rag-h\">Summary</div><div class=\"rag-p\">' + row.rag_result_html + '</div>';
      }

      html += '<div class=\"rag-close\" role=\"button\" tabindex=\"0\">Close</div>';

      pop.innerHTML = html;
      document.body.appendChild(pop);

      var rect = anchorEl.getBoundingClientRect();
      var popRect = pop.getBoundingClientRect();

      var top = window.scrollY + rect.top + 18;
      var left = window.scrollX + rect.left + 18;

      var maxLeft = window.scrollX + document.documentElement.clientWidth - popRect.width - 12;
      var maxTop  = window.scrollY + document.documentElement.clientHeight - popRect.height - 12;

      left = clamp(left, window.scrollX + 12, maxLeft);
      top  = clamp(top,  window.scrollY + 12, maxTop);

      pop.style.left = left + 'px';
      pop.style.top  = top + 'px';

      popPid = pid;

      var closeBtn = pop.querySelector('.rag-close');
      if (closeBtn) {
        closeBtn.addEventListener('click', function(e) { e.stopPropagation(); closePop(); });
        closeBtn.addEventListener('keydown', function(e) {
          if (e.key === 'Enter' || e.key === ' ') { e.preventDefault(); closePop(); }
        });
      }
    }

    function installHandlers() {
      $(document).off('click.rag').on('click.rag', '.rag-info-icon', function(e) {
        e.stopPropagation();
        var pid = $(this).data('pid');
        openPop(this, pid);
      });

      $(document).off('click.rag_out').on('click.rag_out', function() {
        closePop();
      });

      $(window).off('scroll.rag resize.rag').on('scroll.rag resize.rag', function() {
        closePop();
      });
    }

    installStyles();
    installHandlers();
  })();
", rag_payload))

dt_small2 <- datatable(
  df_small2,
  rownames = TRUE,
  escape = FALSE,
  options = list(
    scrollX = FALSE,
    scroller = TRUE,
    pageLength = 25,
    lengthChange = FALSE,
    deferRender = TRUE,
    autoWidth = FALSE,
    columnDefs = list(
      # Visible: 0 rownames, 1 Panel name, 2 Gene count, 4 Info
      # Hidden but searchable: 3 panel_id, 5+ all extra text
      list(visible = FALSE, targets = 3, searchable = TRUE),
      list(visible = FALSE, targets = seq(5, ncol(df_small2)), searchable = TRUE),
      
      list(width = "20px",  targets = 0, orderable = TRUE),
      list(width = "140px", targets = 1),
      list(width = "50px",  targets = 2),
      list(width = "24px",  targets = 4, orderable = FALSE)
    ),
    language = list(
      search = "",
      searchPlaceholder = "Enter natural language query..."
    ),
    initComplete = JS(sprintf("
      function(settings, json) {
        var api = this.api();

        var filter = $('div.dataTables_filter');
        filter.css({'width': '100%%'});

        var input = filter.find('input');
        input.css({
          'width': '100%%',
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

        $('html, body').css({ 'background': 'transparent' });
        $('.dataTables_wrapper').css({ 'background': 'transparent' });

        function normalise(s) {
          if (!s) return '';
          return String(s)
            .toLowerCase()
            .replace(/[\\u2018\\u2019\\u201C\\u201D]/g, ' ')
            .replace(/[,.;:()\\[\\]{}<>\"'`~!@#$%%^&*+=\\\\|\\/\\?-]+/g, ' ')
            .replace(/\\s+/g, ' ')
            .trim();
        }

        function tokenise(q) {
          q = normalise(q);
          if (!q) return [];
          var toks = q.split(' ').filter(Boolean);
          var stop = { and:1, or:1, the:1, a:1, an:1, of:1, to:1, in:1, on:1, for:1, with:1 };
          toks = toks.filter(function(t){ return t.length > 1 && !stop[t]; });

          var uniq = [];
          var seen = {};
          for (var i = 0; i < toks.length; i++) {
            if (!seen[toks[i]]) { seen[toks[i]] = 1; uniq.push(toks[i]); }
          }
          return uniq;
        }

        // cache: dataIndex -> normalised text for that row
        var rowCache = Object.create(null);

        function getRowText(dataIndex, data) {
          if (rowCache[dataIndex]) return rowCache[dataIndex];
          var txt = normalise(data.join(' '));
          rowCache[dataIndex] = txt;
          return txt;
        }

        // install filter once (shared across tables on the page)
        if (!$.fn.dataTable.ext.search._panelapp_nlq_installed) {
          $.fn.dataTable.ext.search.push(function(dtSettings, data, dataIndex) {
            var toks = dtSettings._panelapp_nlq_tokens || [];
            if (toks.length === 0) return true;

            // per-table cache (store on dtSettings so multiple tables do not collide)
            if (!dtSettings._panelapp_rowCache) dtSettings._panelapp_rowCache = Object.create(null);
            var cache = dtSettings._panelapp_rowCache;

            var rowText;
            if (cache[dataIndex]) {
              rowText = cache[dataIndex];
            } else {
              rowText = normalise(data.join(' '));
              cache[dataIndex] = rowText;
            }

            for (var i = 0; i < toks.length; i++) {
              if (rowText.indexOf(toks[i]) === -1) return false;
            }
            return true;
          });
          $.fn.dataTable.ext.search._panelapp_nlq_installed = true;
        }

        // remove DataTables default handler and use debounced draw
        input.off('.DT').off('.panelapp');

        var tmr = null;
        function applyQueryDebounced(raw) {
          clearTimeout(tmr);
          tmr = setTimeout(function() {
            settings._panelapp_nlq_tokens = tokenise(raw);

            // clear cache if you want results to reflect any dynamic changes in row text
            // (normally not needed, but safe)
            settings._panelapp_rowCache = Object.create(null);

            api.draw();
          }, 180);
        }

        input.on('input.panelapp keyup.panelapp search.panelapp', function() {
          applyQueryDebounced(this.value);
        });

        // initial draw state
        settings._panelapp_nlq_tokens = tokenise(input.val());
        api.draw();

        %s
      }
    ", as.character(popover_js)))
  )
) %>%
  formatStyle(
    "Gene count",
    background = styleColorBar(as.numeric(df_small2$`Gene count`), "#43b4eb"),
    backgroundSize = "100% 90%",
    backgroundRepeat = "no-repeat",
    backgroundPosition = "center"
  )

landing_page <- file.path(panelapp_rds_root, "landing_page_beta.html")
saveWidget(dt_small2, landing_page, selfcontained = FALSE)

cat("Saved beta landing page (with AI info icon popover):\n  ", landing_page, "\n")

# The soft match is possible:
# return true if at least ONE token matches, rather than requiring ALL tokens
# However filtering is extremely difficult dues to the large kb on every panel causing many overlaps.
# An option is weighted tokens where not all tokens should be equal. You can give higher weight to:
# gene symbols (all caps + numbers like SERPING1, F12)
# longer words (length ≥ 6)
# rare tokens (only if you want to compute rarity once)
# Then require a minimum score rather than a minimum count. This helps a lot with “large kb overlaps” because common words contribute almost nothing.
# Further, if the user types a multi-word phrase, you can give a bonus if the normalised row text contains the normalised query substring. It helps when users paste text.