# test api ----

# openai_api_test.R

suppressPackageStartupMessages({
  library(httr)
  library(jsonlite)
  library(stringr)
})

api_key <- readLines(
  "~/so/data/keys/switzerlandomics_openai_api_key_panelapprex.txt",
  warn = FALSE
) |>
  str_trim()

if (!nzchar(api_key)) stop("OpenAI API key is empty.")

openai_model <- "gpt-4.1-mini"
# openai_model <- "gpt-5-mini"

body <- list(
  model = openai_model,
  temperature = 0,
  messages = list(
    list(role = "user", content = "Reply with exactly: this is an api test")
  )
)

res <- httr::POST(
  url = "https://api.openai.com/v1/chat/completions",
  httr::add_headers(
    `Content-Type` = "application/json",
    Authorization = paste("Bearer", api_key)
  ),
  body = jsonlite::toJSON(body, auto_unbox = TRUE),
  encode = "raw"
)

status <- httr::status_code(res)
raw_text <- httr::content(res, as = "text", encoding = "UTF-8")

cat("HTTP status:", status, "\n")

if (status != 200L) {
  err <- tryCatch(jsonlite::fromJSON(raw_text), error = function(e) NULL)
  if (!is.null(err$error)) {
    cat("Error code:", err$error$code %||% "NA", "\n")
    cat("Message:", err$error$message %||% "NA", "\n")
  } else {
    cat("Raw response:\n", raw_text, "\n")
  }
  stop("OpenAI request failed.")
}

parsed <- jsonlite::fromJSON(raw_text)
cat("Model:", parsed$model, "\n")
cat("Response:\n")
cat(parsed$choices$message$content, "\n")

# cat(parsed$choices[[0]]$message$content, "\n")

# `%||%` <- function(x, y) if (is.null(x) || is.na(x) || !nzchar(as.character(x))) y else x

# end test ----

# panel_rag_uniprot_gene_level.R

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(purrr)
  library(stringr)
  library(httr)
  library(jsonlite)
})
  
# .................................................
# Config ----
# .................................................
panelapp_rds_root <- "../data"
panelapprex_root  <- "~/mnt/atlas_data_big/data/panelapprex"
# PanelApp combined data (gene level, multiple panels)
path_panels_rds <- file.path(panelapp_rds_root, "PanelAppData_genes_combined_Rds")

uniprot_dir <- file.path(panelapprex_root, "uniprot_up000005640")
uniprot_dat <- file.path(uniprot_dir, "UP000005640_9606.dat.gz")

out_dir_processed <- "~/mnt/atlas_data_big/data/panelapprex/uniprot_kb_up000005640_processed"

# Output directory for RAG panel summaries
out_dir_rag <- file.path(panelapprex_root, "uniprot_kb_up000005640_processed_rag")
dir.create(out_dir_rag, showWarnings = FALSE, recursive = TRUE)

openai_model <- "gpt-4.1-mini" # 400,000 TPM / 200 RPM / 4,000,000 TPD
# openai_model <- "gpt-5-nano"
# openai_model <- "gpt-5.1-mini-2025-08-07"
# openai_model <- "gpt-4.1" # smartest non-reasoning model

# Switch between single test panel and all panels
process_all_panels <- FALSE       # set TRUE later for full run
process_all_panels <- TRUE       # set TRUE later for full run
# test_panel_id      <- 398L        # used only if process_all_panels == FALSE
# test_panel_id      <- 285L        # size > 2000
# test_panel_id      <- 467L        # size ~ 1000
# test_panel_id      <- 3L        # size < 10

# per gene truncation and context size
skip_if_genes_gt <- 1000L
max_chars_per_gene <- 99999L # remove for full run
max_genes_per_panel <- 9999L # remove for full run

# OpenAI key and model
api_key <- readLines(
  "~/so/data/keys/switzerlandomics_openai_api_key_panelapprex.txt",
  warn = FALSE
) |>
  str_trim()

if (!nzchar(api_key)) {
  stop("OpenAI API key is empty. Check the key file path.")
}


# PanelApp combined data (gene level, multiple panels)
# path_panels_rds <- file.path(path_data, "PanelAppData_genes_combined_Rds")

# UniProtKB gene level table from build_uniprotkb_up000005640_whole_proteome.R
path_uniprot_gene_rds <- file.path(
  # path_data,
  # uniprot_dir,
  out_dir_processed,
  # "uniprot_kb_up000005640",
  "uniprot_gene_level.rds"
)


# .................................................
# Load data ----
# .................................................

if (!file.exists(path_panels_rds)) {
  stop("PanelApp RDS not found at: ", path_panels_rds)
}

if (!file.exists(path_uniprot_gene_rds)) {
  stop("UniProt gene level RDS not found at: ", path_uniprot_gene_rds)
}

df_panels_raw <- readRDS(path_panels_rds)
df_uniprot_raw <- readRDS(path_uniprot_gene_rds)

# force to tibbles so dplyr verbs behave predictably
df_panels <- as_tibble(df_panels_raw)
df_uniprot_gene <- as_tibble(df_uniprot_raw)

# Basic sanity checks
required_panel_cols <- c(
  "panel_id",
  "entity_name",
  "name",
  "disease_group",
  "disease_sub_group",
  "mode_of_inheritance",
  "phenotypes"
)

missing_panel <- setdiff(required_panel_cols, names(df_panels))
if (length(missing_panel) > 0) {
  stop("df_panels is missing columns: ", paste(missing_panel, collapse = ", "))
}

required_uniprot_cols <- c("SYMBOL", "tier1_text")
missing_uniprot <- setdiff(required_uniprot_cols, names(df_uniprot_gene))
if (length(missing_uniprot) > 0) {
  stop("df_uniprot_gene is missing columns: ", paste(missing_uniprot, collapse = ", "))
}

# Decide which panels to process
if (process_all_panels) {
  panel_ids <- sort(unique(df_panels$panel_id))
} else {
  panel_ids <- test_panel_id
}

cat("Panels to process:", paste(panel_ids, collapse = ", "), "\n")

# .................................................
# reduce content for API ----
# .................................................

clean_uniprot_for_rag <- function(txt) {
  txt <- ifelse(is.na(txt), "", txt)
  
  # Remove SIMILARITY
  txt <- stringr::str_remove_all(
    txt,
    "SIMILARITY:\\s+.*?(?=\\|\\||$)"
  )
  
  # # Remove SUBCELLULAR LOCATION
  # txt <- stringr::str_remove_all(
  #   txt,
  #   "SUBCELLULAR LOCATION:\\s+.*?(?=\\|\\||$)"
  # )
  
  txt <- stringr::str_remove_all(
    txt,
    "Note=The disease is caused by variants affecting the gene represented in this entry\\."
  )
  
  # Remove ECO tags
  txt <- stringr::str_remove_all(
    txt,
    "\\{ECO:[^}]+\\}"
  )
  
  # Remove PubMed references
  txt <- stringr::str_remove_all(
    txt,
    "PubMed:\\d+"
  )
  
  # Remove "(By similarity)"
  txt <- stringr::str_remove_all(
    txt,
    "\\(By similarity\\)"
  )
  
  # Remove footer
  txt <- stringr::str_remove_all(
    txt,
    "Copyrighted by the UniProt Consortium.*$"
  )
  
  txt <- stringr::str_squish(txt)
  txt
}


# .................................................
# Static prompts
# .................................................

system_prompt <- paste(
  "Your task is to compress mechanistic biological knowledge into a brief, high-signal summary",
  "that enables clinical and mechanistic search.",
  "Output must be deterministic, concise, and medically meaningful.",
  "Use language that matches how clinicians describe patients, including key syndromes and phenotype phrases.",
  "Do not include filler words, disclaimers, or explanations.",
  "Use 3 to 6 bullet points.",
  "Each bullet must express a mechanistic pathway, disease mechanism, typical phenotypes, or organ involvement.",
  "Do not mention gene names.",
  "Do not mention the term 'panel' or describe the task.",
  "Do not repeat yourself.",
  "Do not invent facts.",
  "Your goal is to maximise the semantic search value of the output for clinicians.",
  sep = " "
)

# .................................................
# Main loop over panels
# .................................................


for (pid in panel_ids) {
  cat("\n=== Processing panel", pid, "===\n")
  
  rds_path <- file.path(out_dir_rag, paste0("panel_", pid, "_rag_summary.rds"))
  
  if (file.exists(rds_path)) {
    cat("Existing summary found for panel", pid, "skipping.\n")
    next
  }
  
  txt_full_path <- file.path(out_dir_rag, paste0("panel_", pid, "_rag_full.txt"))
  txt_res_path <- file.path(out_dir_rag, paste0("panel_", pid, "_rag_result.txt"))
  txt_over_path <- file.path(out_dir_rag, paste0("panel_", pid, "_rag_overview.txt"))
  
  # rest of the loop as before, but remove the second rds_path / txt_* redefinition
# for (pid in panel_ids) {
#   cat("\n=== Processing panel", pid, "===\n")
#   
#   rds_path <- file.path(out_dir_rag, paste0("panel_", pid, "_rag_summary.rds"))
#   txt_path <- file.path(out_dir_rag, paste0("panel_", pid, "_rag_summary.txt"))
#   
#   # Skip if already done
#   if (file.exists(rds_path) && file.exists(txt_path)) {
#     cat("Existing summary found for panel", pid, "- skipping.\n")
#     next
#   }
  
  # Extract gene set for this panel
  panel_genes <- df_panels |>
    filter(panel_id == pid) |>
    distinct(entity_name, .keep_all = TRUE) |>
    transmute(
      panel_id,
      panel_name = name,
      disease_group,
      disease_sub_group,
      mode_of_inheritance,
      phenotypes,
      SYMBOL = entity_name
    )
  
  if (nrow(panel_genes) == 0) {
    cat("No records found for panel_id =", pid, "- skipping.\n")
    next
  }
  
  # skip massive panels
  n_panel_genes <- nrow(panel_genes)
  if (n_panel_genes > skip_if_genes_gt) {
    cat("Skipping panel", pid, "because gene count",
        n_panel_genes, "exceeds", skip_if_genes_gt, "\n")
    next
  }
  
  panel_join <- panel_genes |>
    left_join(df_uniprot_gene, by = "SYMBOL")
  
  # Keep only genes with non empty tier1_text and truncate for token safety
  panel_join <- panel_join |>
    mutate(
      # tier1_text = ifelse(is.na(tier1_text), "", tier1_text),
      # tier1_text = str_squish(tier1_text),
      # tier1_text = ifelse(
      #   nchar(tier1_text) > max_chars_per_gene,
      #   paste0(substr(tier1_text, 1L, max_chars_per_gene), " ..."),
      #   tier1_text
      
        tier1_text = ifelse(is.na(tier1_text), "", tier1_text),
        tier1_text = clean_uniprot_for_rag(tier1_text),
        tier1_text = str_squish(tier1_text),
        tier1_text = ifelse(
          nchar(tier1_text) > max_chars_per_gene,
          paste0(substr(tier1_text, 1L, max_chars_per_gene), " ..."),
          tier1_text
        )
    ) |>
    filter(nzchar(tier1_text))
  
  if (nrow(panel_join) == 0) {
    cat("No UniProt tier1_text available for panel", pid, "- skipping.\n")
    next
  }
  
  # Optionally cap number of genes for context (use base indexing to avoid slice/Rle issues)
  # n_keep <- min(nrow(panel_join), max_genes_per_panel)
  # panel_join_small <- panel_join[seq_len(n_keep), , drop = FALSE]
  n_keep <- min(nrow(panel_join), max_genes_per_panel, skip_if_genes_gt)
  panel_join_small <- panel_join[seq_len(n_keep), , drop = FALSE]
  
  # Build gene context text
  gene_context_text <- panel_join_small |>
    mutate(
      phenotype_snip = ifelse(
        !is.na(phenotypes) & nzchar(phenotypes),
        paste0("Reported phenotypes: ", phenotypes, ". "),
        ""
      ),
      snippet = paste0(
        "Gene ", SYMBOL, ": ",
        phenotype_snip,
        tier1_text
      )
    ) |>
    pull(snippet) |>
    paste(collapse = "\n")
  
  # Panel metadata text
  panel_meta <- panel_join_small |>
    distinct(
      panel_id,
      panel_name,
      disease_group,
      disease_sub_group,
      mode_of_inheritance
    )
  
  panel_meta <- as.data.frame(panel_meta)[1, , drop = FALSE]
  
  meta_text <- panel_meta |>
    transmute(
      txt = paste0(
        "Panel id: ", panel_id, "\n",
        "Panel name: ", panel_name, "\n",
        "Disease group: ", disease_group, "\n",
        "Disease subgroup: ", disease_sub_group, "\n",
        "Modes of inheritance in this panel: ", mode_of_inheritance, "\n"
      )
    ) |>
    pull(txt)
  
  # Build user prompt for this panel
  user_prompt <- paste0(
    "You will receive mechanistic summaries for all genes in one diagnostic gene set.\n\n",
    "Follow these formatting rules exactly:\n\n",
    "1. First output the line 'PANELAPPREX_AI_ANALYSIS'.\n",
    "2. On the next line, output 'GENES: ' followed by a comma separated list of all gene symbols you considered (for example: GENES: ABL1, ACTA2, SMAD3).\n",
    "3. Then write one or two short paragraphs explaining which shared biological mechanisms, pathways, disease processes, and characteristic phenotypes you identified from the gene summaries.\n",
    "4. After the analysis, output the line 'PANELAPPREX_AI_OVERVIEW' on its own line.\n",
    "5. Under PANELAPPREX_AI_OVERVIEW, write a single concise paragraph that a clinician could read to understand the main clinical and phenotypic focus of this gene set. Use phenotype and syndrome language, not gene symbols.\n",
    "6. After the overview, output the line 'PANELAPPREX_AI_RESULT' on its own line.\n",
    "7. Under PANELAPPREX_AI_RESULT, output 3 to 6 bullet points with the final condensed summary.\n",
    "   - Focus on mechanistic signals and characteristic phenotypes.\n",
    "   - Do not mention gene names in these RESULT bullets.\n",
    "   - Do not include any other text before or after these three sections.\n\n",
    "Here is the metadata:\n\n",
    meta_text,
    "\n",
    "Here is the input text from all genes in the gene set:\n\n",
    gene_context_text
  )
  
  body <- list(
    model = openai_model,
    temperature = 0.1,
    messages = list(
      list(
        role = "system",
        content = system_prompt
      ),
      list(
        role = "user",
        content = user_prompt
      )
    )
  )
  
  # Call OpenAI with basic error handling so a single failure does not kill the batch
  # ok <- TRUE
  # summary_text <- NULL
  # 
  # tryCatch(
  #   {
  #     res <- POST(
  #       url = "https://api.openai.com/v1/chat/completions",
  #       add_headers(
  #         `Content-Type` = "application/json",
  #         Authorization = paste("Bearer", api_key)
  #       ),
  #       body = toJSON(body, auto_unbox = TRUE)
  #     )
  #     
  #     stop_for_status(res)
  #     
  #     parsed <- content(res, as = "parsed", encoding = "UTF-8")
  #     
  #     if (is.null(parsed$choices) || length(parsed$choices) < 1) {
  #       stop("OpenAI response did not contain any choices")
  #     }
  #     
  #     summary_text <- parsed$choices[[1]]$message$content
  #   },
  #   error = function(e) {
  #     ok <<- FALSE
  #     cat("OpenAI call failed for panel", pid, ":", conditionMessage(e), "\n")
  #   }
  # )
  # Call OpenAI with basic error handling so a single failure does not kill the batch
  ok <- TRUE
  summary_text <- NULL
  
  tryCatch(
    {
      res <- POST(
        url = "https://api.openai.com/v1/chat/completions",
        add_headers(
          `Content-Type` = "application/json",
          Authorization = paste("Bearer", api_key)
        ),
        body = toJSON(body, auto_unbox = TRUE)
      )
      
      status <- httr::status_code(res)
      if (status != 200L) {
        raw_text <- httr::content(res, as = "text", encoding = "UTF-8")
        err_obj <- NULL
        err_code <- NA_character_
        err_msg  <- NA_character_
        
        err_obj <- tryCatch(
          jsonlite::fromJSON(raw_text),
          error = function(e) NULL
        )
        
        if (!is.null(err_obj$error$code)) {
          err_code <- err_obj$error$code
        }
        if (!is.null(err_obj$error$message)) {
          err_msg <- err_obj$error$message
        }
        
        cat(
          "OpenAI HTTP error for panel", pid,
          "- status:", status,
          "code:", ifelse(is.na(err_code), "NA", err_code), "\n"
        )
        if (!is.na(err_msg)) {
          cat("Message:", err_msg, "\n")
        } else {
          cat("Raw response body:\n", raw_text, "\n")
        }
        
        httr::stop_for_status(res)
      }
      
      parsed <- httr::content(res, as = "parsed", encoding = "UTF-8")
      
      if (is.null(parsed$choices) || length(parsed$choices) < 1) {
        stop("OpenAI response did not contain any choices")
      }
      
      summary_text <- parsed$choices[[1]]$message$content
    },
    error = function(e) {
      ok <<- FALSE
      cat("OpenAI call failed for panel", pid, ":", conditionMessage(e), "\n")
    }
  )
  
  if (!ok || is.null(summary_text)) {
    cat("Skipping save for panel", pid, "due to OpenAI error.\n")
    next
  }
  
  # Split into analysis, overview, and result using the delimiters
  parts1 <- stringr::str_split_fixed(summary_text, "PANELAPPREX_AI_OVERVIEW", n = 2)
  
  analysis_text <- stringr::str_trim(parts1[, 1])
  tail_text     <- if (ncol(parts1) >= 2) stringr::str_trim(parts1[, 2]) else ""
  
  parts2 <- stringr::str_split_fixed(tail_text, "PANELAPPREX_AI_RESULT", n = 2)
  
  overview_text <- stringr::str_trim(parts2[, 1])
  result_text   <- if (ncol(parts2) >= 2) stringr::str_trim(parts2[, 2]) else ""
  
  # Safety checks
  if (!stringr::str_detect(analysis_text, "PANELAPPREX_AI_ANALYSIS")) {
    warning("Response for panel ", pid, " does not contain PANELAPPREX_AI_ANALYSIS as expected.")
  }
  if (!nzchar(overview_text)) {
    warning("Response for panel ", pid, " does not contain a non empty PANELAPPREX_AI_OVERVIEW section.")
  }
  if (!nzchar(result_text)) {
    warning("Response for panel ", pid, " does not contain a non empty PANELAPPREX_AI_RESULT section.")
  }
  
  # Save outputs
  summary_obj <- list(
    panel_id = panel_meta$panel_id[[1]],
    panel_name = panel_meta$panel_name[[1]],
    disease_group = panel_meta$disease_group[[1]],
    disease_sub_group = panel_meta$disease_sub_group[[1]],
    mode_of_inheritance = panel_meta$mode_of_inheritance[[1]],
    analysis_text = analysis_text,
    overview_text = overview_text,
    summary_text = result_text,
    raw_text = summary_text,
    model = openai_model,
    n_genes_context = nrow(panel_join_small),
    max_chars_per_gene = max_chars_per_gene
  )
  
  # rds_path        <- file.path(out_dir_rag, paste0("panel_", pid, "_rag_summary.rds"))
  # txt_full_path   <- file.path(out_dir_rag, paste0("panel_", pid, "_rag_full.txt"))
  # txt_res_path    <- file.path(out_dir_rag, paste0("panel_", pid, "_rag_result.txt"))
  # txt_over_path   <- file.path(out_dir_rag, paste0("panel_", pid, "_rag_overview.txt"))
  
  saveRDS(summary_obj, rds_path)
  writeLines(summary_text,  txt_full_path)   # full reasoning + overview + result
  writeLines(overview_text, txt_over_path)   # overview paragraph only
  writeLines(result_text,   txt_res_path)    # final bullets only
  
  cat(
    "Saved RAG summary for panel",
    pid,
    "to:\n  ", rds_path,
    "\n  ", txt_full_path,
    "\n  ", txt_over_path,
    "\n  ", txt_res_path, "\n"
  )
}

cat("Character length of prompt:", nchar(user_prompt), "\n")
summary_text
# head(user_prompt)
cat("\nDone.\n")
