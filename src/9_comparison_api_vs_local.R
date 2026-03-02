#!/usr/bin/env Rscript

set.seed(666)

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(stringr)
  library(jsonlite)
  library(httr)
  library(readr)
  library(tibble)
})

# ------------------------------------------------
# Config
# ------------------------------------------------

seed_panels <- 398L
n_random_extra_panels <- 149L
random_seed <- 1L
case_index <- 1L

base_url <- "https://panelapp.genomicsengland.co.uk/api/v1"

path_data <- "../data"
dir.create(path_data, showWarnings = FALSE, recursive = TRUE)

path_benchmark <- file.path(path_data, "benchmark.tsv")
path_core_rds  <- file.path(path_data, "PanelAppData_genes_combined_Rds")
path_api_key   <- "../credentials/creds_api_key.txt"

results_tag <- sprintf(
  "cold_start_seed%s_plus%d_total%s_case%s",
  paste(seed_panels, collapse = "-"),
  n_random_extra_panels,
  "%s",
  case_index
)

path_results_tsv <- file.path(path_data, sprintf("benchmark_results_%s.tsv", results_tag))
path_subset_tsv  <- file.path(path_data, sprintf("benchmark_subset_%s.tsv", results_tag))
path_summary_tsv <- file.path(path_data, sprintf("benchmark_summary_%s.tsv", results_tag))
path_problems_tsv <- file.path(path_data, sprintf("benchmark_ge_problems_%s.tsv", results_tag))

# ------------------------------------------------
# Helpers
# ------------------------------------------------

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0) return(y)
  x
}

stop_if_missing <- function(path) {
  if (!file.exists(path)) stop("Missing file: ", path, call. = FALSE)
}

read_api_key <- function(path) {
  if (!file.exists(path)) return(NA_character_)
  x <- readLines(path, warn = FALSE)
  x <- trimws(x)
  x <- x[nzchar(x)]
  if (length(x) == 0) return(NA_character_)
  x[[1]]
}

normalise <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  x <- tolower(x)
  x <- str_replace_all(x, "[\u2018\u2019\u201C\u201D]", " ")
  x <- str_replace_all(x, "[,.;:()\\[\\]{}<>\"'`~!@#$%^&*+=\\\\|/\\?-]+", " ")
  x <- str_replace_all(x, "\\s+", " ")
  x <- str_trim(x)
  x
}

tokenise <- function(query) {
  x <- normalise(query)
  if (length(x) == 0 || !nzchar(x[[1]])) return(character())
  toks <- strsplit(x[[1]], " ", fixed = TRUE)[[1]]
  toks <- toks[nzchar(toks)]
  stop_words <- c("and", "or", "the", "a", "an", "of", "to", "in", "on", "for", "with")
  toks <- toks[nchar(toks) > 1 & !(toks %in% stop_words)]
  unique(toks)
}

safe_text <- function(x) {
  if (is.null(x) || length(x) == 0) return("")
  x <- unlist(x, use.names = FALSE)
  x <- as.character(x)
  x <- x[!is.na(x)]
  x <- trimws(x)
  x <- x[nzchar(x)]
  if (length(x) == 0) return("")
  paste(unique(x), collapse = "; ")
}

first_nonempty <- function(x) {
  x <- as.character(x)
  x <- x[!is.na(x)]
  x <- trimws(x)
  x <- x[nzchar(x)]
  if (length(x) == 0) return(NA_character_)
  x[[1]]
}

coerce_col <- function(df, col) {
  if (!col %in% names(df)) return(rep("", nrow(df)))
  vapply(df[[col]], safe_text, character(1))
}

extract_case_col <- function(df) {
  hit <- names(df)[tolower(names(df)) %in% c("case_study", "case study")]
  if (length(hit) == 0) stop("Could not find case study column in benchmark.tsv", call. = FALSE)
  hit[[1]]
}

extract_query_col <- function(df) {
  hit <- names(df)[tolower(names(df)) %in% c("query", "query_terms")]
  if (length(hit) == 0) stop("Could not find query column in benchmark.tsv", call. = FALSE)
  hit[[1]]
}

extract_causal_col <- function(df) {
  hit <- names(df)[tolower(names(df)) %in% c("causal_gene", "causal gene")]
  if (length(hit) == 0) stop("Could not find causal gene column in benchmark.tsv", call. = FALSE)
  hit[[1]]
}

score_panel_text <- function(search_text_norm, tokens, causal_gene, entity_names) {
  token_score <- if (length(tokens) == 0) 0L else {
    sum(vapply(tokens, function(tok) str_detect(search_text_norm, fixed(tok)), logical(1)))
  }
  
  causal_hit <- as.integer(str_detect(
    normalise(entity_names %||% ""),
    fixed(normalise(causal_gene))
  ))
  
  list(
    token_score = as.integer(token_score),
    causal_gene_in_panel = causal_hit
  )
}

rank_panels <- function(index_df, query, causal_gene) {
  tokens <- tokenise(query)
  
  if (nrow(index_df) == 0) {
    return(list(
      ranked = index_df,
      tokens = tokens
    ))
  }
  
  scored <- index_df |>
    rowwise() |>
    mutate(
      token_score = score_panel_text(search_text_norm, tokens, causal_gene, entity_names)$token_score,
      causal_gene_in_panel = score_panel_text(search_text_norm, tokens, causal_gene, entity_names)$causal_gene_in_panel
    ) |>
    ungroup() |>
    arrange(desc(token_score), desc(causal_gene_in_panel), gene_count, panel_id)
  
  list(
    ranked = scored,
    tokens = tokens
  )
}

summarise_search_result <- function(
    ranked_df,
    query,
    causal_gene,
    method,
    subset_size,
    elapsed_total_seconds,
    n_api_calls,
    n_panels_requested,
    n_panels_indexed
) {
  tokens <- tokenise(query)
  
  if (nrow(ranked_df) == 0) {
    return(tibble(
      method = method,
      subset_size = subset_size,
      query = query,
      query_tokens = paste(tokens, collapse = " "),
      n_query_tokens = length(tokens),
      elapsed_total_seconds = elapsed_total_seconds,
      top_panel_id = NA_integer_,
      top_panel_name = NA_character_,
      n_candidate_panels_returned = 0L,
      causal_gene = causal_gene,
      causal_gene_found_in_any_returned_panel = 0L,
      best_panel_contains_causal_gene = 0L,
      n_api_calls = as.integer(n_api_calls),
      n_panels_requested = as.integer(n_panels_requested),
      n_panels_indexed = as.integer(n_panels_indexed)
    ))
  }
  
  returned <- ranked_df |>
    filter(token_score > 0)
  
  top <- ranked_df[1, , drop = FALSE]
  
  tibble(
    method = method,
    subset_size = subset_size,
    query = query,
    query_tokens = paste(tokens, collapse = " "),
    n_query_tokens = length(tokens),
    elapsed_total_seconds = elapsed_total_seconds,
    top_panel_id = as.integer(top$panel_id[[1]]),
    top_panel_name = top$panel_name[[1]],
    n_candidate_panels_returned = nrow(returned),
    causal_gene = causal_gene,
    causal_gene_found_in_any_returned_panel = as.integer(any(returned$causal_gene_in_panel == 1L)),
    best_panel_contains_causal_gene = as.integer(top$causal_gene_in_panel[[1]]),
    n_api_calls = as.integer(n_api_calls),
    n_panels_requested = as.integer(n_panels_requested),
    n_panels_indexed = as.integer(n_panels_indexed)
  )
}

empty_ge_index <- function() {
  tibble(
    panel_id = integer(),
    panel_name = character(),
    disease_group = character(),
    disease_sub_group = character(),
    gene_count = integer(),
    entity_names = character(),
    phenotypes = character(),
    mode_of_inheritance = character(),
    mode_of_pathogenicity = character(),
    publications = character(),
    evidence = character(),
    tags = character(),
    alias = character(),
    alias_name = character(),
    gene_name = character(),
    gene_symbol = character(),
    hgnc_symbol = character(),
    omim_gene = character(),
    hgnc_id = character(),
    ensembl_location = character(),
    ensembl_id = character(),
    search_text = character(),
    search_text_norm = character(),
    token_score = integer(),
    causal_gene_in_panel = integer()
  ) |>
    select(-token_score, -causal_gene_in_panel)
}

api_get_json <- function(url, api_key = NA_character_) {
  headers <- httr::add_headers(Accept = "application/json")
  if (!is.na(api_key) && nzchar(api_key)) {
    headers <- httr::add_headers(
      Accept = "application/json",
      `X-CSRFToken` = api_key
    )
  }
  
  res <- httr::GET(url, headers)
  status <- httr::status_code(res)
  
  if (status >= 200L && status < 300L) {
    txt <- httr::content(res, as = "text", encoding = "UTF-8")
    return(list(
      ok = TRUE,
      status_code = status,
      obj = jsonlite::fromJSON(txt, flatten = TRUE),
      error_message = NA_character_
    ))
  }
  
  txt <- tryCatch(
    httr::content(res, as = "text", encoding = "UTF-8"),
    error = function(e) ""
  )
  
  list(
    ok = FALSE,
    status_code = status,
    obj = NULL,
    error_message = paste0("HTTP ", status, if (nzchar(txt)) paste0(": ", txt) else "")
  )
}

fetch_all_ge_panel_metadata <- function(api_key = NA_character_) {
  page <- 1L
  out <- list()
  total_api_calls <- 0L
  
  repeat {
    got <- api_get_json(sprintf("%s/panels/?page=%d", base_url, page), api_key = api_key)
    total_api_calls <- total_api_calls + 1L
    
    if (!isTRUE(got$ok)) {
      stop("Failed to fetch GE metadata page ", page, ". ", got$error_message, call. = FALSE)
    }
    
    obj <- got$obj
    res <- obj$results %||% NULL
    if (is.null(res) || nrow(as.data.frame(res)) == 0) break
    
    out[[length(out) + 1L]] <- as_tibble(as.data.frame(res, stringsAsFactors = FALSE))
    
    next_url <- obj$`next` %||% NA_character_
    if (is.na(next_url) || !nzchar(next_url)) break
    
    page <- page + 1L
  }
  
  list(
    df = bind_rows(out),
    api_calls = total_api_calls
  )
}

choose_panel_subset <- function(meta_df, seed_panels, n_random_extra_panels = 0L, random_seed = 1L) {
  all_ids <- sort(unique(as.integer(meta_df$id)))
  all_ids <- all_ids[!is.na(all_ids)]
  
  seed_panels <- unique(as.integer(seed_panels))
  seed_panels <- seed_panels[seed_panels %in% all_ids]
  
  if (length(seed_panels) == 0) {
    stop("None of the seed panels were found in GE panel metadata.", call. = FALSE)
  }
  
  pool <- setdiff(all_ids, seed_panels)
  
  set.seed(random_seed)
  extra <- if (n_random_extra_panels > 0L) {
    sample(pool, size = min(n_random_extra_panels, length(pool)), replace = FALSE)
  } else {
    integer()
  }
  
  sort(unique(c(seed_panels, extra)))
}

collapse_panelapp_subset <- function(df_core_subset) {
  df_select <- tibble(
    panel_id = as.integer(df_core_subset$panel_id),
    panel_name = as.character(df_core_subset$name %||% ""),
    entity_name = as.character(df_core_subset$entity_name %||% ""),
    phenotypes = coerce_col(df_core_subset, "phenotypes"),
    mode_of_inheritance = coerce_col(df_core_subset, "mode_of_inheritance"),
    disease_group = coerce_col(df_core_subset, "disease_group"),
    disease_sub_group = coerce_col(df_core_subset, "disease_sub_group"),
    ensembl_location = coerce_col(df_core_subset, "gene_data.ensembl_genes.GRch38.90.location"),
    ensembl_id = coerce_col(df_core_subset, "gene_data.ensembl_genes.GRch38.90.ensembl_id"),
    gene_name = coerce_col(df_core_subset, "gene_data.gene_name"),
    gene_symbol = coerce_col(df_core_subset, "gene_data.gene_symbol"),
    hgnc_symbol = coerce_col(df_core_subset, "gene_data.hgnc_symbol"),
    omim_gene = coerce_col(df_core_subset, "gene_data.omim_gene"),
    hgnc_id = coerce_col(df_core_subset, "gene_data.hgnc_id")
  )
  
  df_select |>
    group_by(panel_id) |>
    summarise(
      panel_name = first_nonempty(panel_name),
      disease_group = paste(unique(disease_group[nzchar(disease_group)]), collapse = "; "),
      disease_sub_group = paste(unique(disease_sub_group[nzchar(disease_sub_group)]), collapse = "; "),
      gene_count = n_distinct(entity_name[nzchar(entity_name)]),
      entity_names = paste(unique(entity_name[nzchar(entity_name)]), collapse = "; "),
      phenotypes = paste(unique(phenotypes[nzchar(phenotypes)]), collapse = "; "),
      mode_of_inheritance = paste(unique(mode_of_inheritance[nzchar(mode_of_inheritance)]), collapse = "; "),
      ensembl_location = paste(unique(ensembl_location[nzchar(ensembl_location)]), collapse = "; "),
      ensembl_id = paste(unique(ensembl_id[nzchar(ensembl_id)]), collapse = "; "),
      gene_name = paste(unique(gene_name[nzchar(gene_name)]), collapse = "; "),
      gene_symbol = paste(unique(gene_symbol[nzchar(gene_symbol)]), collapse = "; "),
      hgnc_symbol = paste(unique(hgnc_symbol[nzchar(hgnc_symbol)]), collapse = "; "),
      omim_gene = paste(unique(omim_gene[nzchar(omim_gene)]), collapse = "; "),
      hgnc_id = paste(unique(hgnc_id[nzchar(hgnc_id)]), collapse = "; "),
      .groups = "drop"
    ) |>
    mutate(
      search_text = paste(
        panel_name,
        disease_group,
        disease_sub_group,
        entity_names,
        phenotypes,
        mode_of_inheritance,
        ensembl_location,
        ensembl_id,
        gene_name,
        gene_symbol,
        hgnc_symbol,
        omim_gene,
        hgnc_id,
        sep = " ; "
      ),
      search_text_norm = normalise(search_text)
    )
}

fetch_ge_panel_detail <- function(panel_id, api_key = NA_character_) {
  got <- api_get_json(sprintf("%s/panels/%s/", base_url, panel_id), api_key = api_key)
  
  if (isTRUE(got$ok)) {
    return(list(
      status = "ok",
      obj = got$obj,
      api_calls = 1L,
      message = NA_character_
    ))
  }
  
  if (identical(got$status_code, 404L)) {
    return(list(
      status = "missing_panel",
      obj = NULL,
      api_calls = 1L,
      message = got$error_message
    ))
  }
  
  list(
    status = "fetch_error",
    obj = NULL,
    api_calls = 1L,
    message = got$error_message
  )
}

build_ge_single_panel_index <- function(detail_obj) {
  genes <- detail_obj$genes %||% NULL
  if (is.null(genes) || !is.data.frame(genes) || nrow(genes) == 0) {
    stop("No genes found in GE panel detail response for panel ", detail_obj$id %||% NA, call. = FALSE)
  }
  
  gene_entity_name <- if ("entity_name" %in% names(genes)) as.character(genes$entity_name) else rep("", nrow(genes))
  gene_phenotypes  <- if ("phenotypes" %in% names(genes)) vapply(genes$phenotypes, safe_text, character(1)) else rep("", nrow(genes))
  gene_moi         <- if ("mode_of_inheritance" %in% names(genes)) as.character(genes$mode_of_inheritance) else rep("", nrow(genes))
  gene_mop         <- if ("mode_of_pathogenicity" %in% names(genes)) as.character(genes$mode_of_pathogenicity) else rep("", nrow(genes))
  gene_pubs        <- if ("publications" %in% names(genes)) vapply(genes$publications, safe_text, character(1)) else rep("", nrow(genes))
  gene_evidence    <- if ("evidence" %in% names(genes)) vapply(genes$evidence, safe_text, character(1)) else rep("", nrow(genes))
  gene_tags        <- if ("tags" %in% names(genes)) vapply(genes$tags, safe_text, character(1)) else rep("", nrow(genes))
  
  alias_col        <- if ("gene_data.alias" %in% names(genes)) vapply(genes[["gene_data.alias"]], safe_text, character(1)) else rep("", nrow(genes))
  alias_name_col   <- if ("gene_data.alias_name" %in% names(genes)) vapply(genes[["gene_data.alias_name"]], safe_text, character(1)) else rep("", nrow(genes))
  gene_name_col    <- if ("gene_data.gene_name" %in% names(genes)) as.character(genes[["gene_data.gene_name"]]) else rep("", nrow(genes))
  gene_symbol_col  <- if ("gene_data.gene_symbol" %in% names(genes)) as.character(genes[["gene_data.gene_symbol"]]) else rep("", nrow(genes))
  hgnc_symbol_col  <- if ("gene_data.hgnc_symbol" %in% names(genes)) as.character(genes[["gene_data.hgnc_symbol"]]) else rep("", nrow(genes))
  omim_gene_col    <- if ("gene_data.omim_gene" %in% names(genes)) vapply(genes[["gene_data.omim_gene"]], safe_text, character(1)) else rep("", nrow(genes))
  hgnc_id_col      <- if ("gene_data.hgnc_id" %in% names(genes)) as.character(genes[["gene_data.hgnc_id"]]) else rep("", nrow(genes))
  grch38_loc_col   <- if ("gene_data.ensembl_genes.GRch38.90.location" %in% names(genes)) as.character(genes[["gene_data.ensembl_genes.GRch38.90.location"]]) else rep("", nrow(genes))
  grch38_id_col    <- if ("gene_data.ensembl_genes.GRch38.90.ensembl_id" %in% names(genes)) as.character(genes[["gene_data.ensembl_genes.GRch38.90.ensembl_id"]]) else rep("", nrow(genes))
  
  tibble(
    panel_id = as.integer(detail_obj$id),
    panel_name = as.character(detail_obj$name %||% ""),
    disease_group = as.character(detail_obj$disease_group %||% ""),
    disease_sub_group = as.character(detail_obj$disease_sub_group %||% ""),
    gene_count = length(unique(gene_entity_name[nzchar(gene_entity_name)])),
    entity_names = paste(unique(gene_entity_name[nzchar(gene_entity_name)]), collapse = "; "),
    phenotypes = paste(unique(gene_phenotypes[nzchar(gene_phenotypes)]), collapse = "; "),
    mode_of_inheritance = paste(unique(gene_moi[nzchar(gene_moi)]), collapse = "; "),
    mode_of_pathogenicity = paste(unique(gene_mop[nzchar(gene_mop)]), collapse = "; "),
    publications = paste(unique(gene_pubs[nzchar(gene_pubs)]), collapse = "; "),
    evidence = paste(unique(gene_evidence[nzchar(gene_evidence)]), collapse = "; "),
    tags = paste(unique(gene_tags[nzchar(gene_tags)]), collapse = "; "),
    alias = paste(unique(alias_col[nzchar(alias_col)]), collapse = "; "),
    alias_name = paste(unique(alias_name_col[nzchar(alias_name_col)]), collapse = "; "),
    gene_name = paste(unique(gene_name_col[nzchar(gene_name_col)]), collapse = "; "),
    gene_symbol = paste(unique(gene_symbol_col[nzchar(gene_symbol_col)]), collapse = "; "),
    hgnc_symbol = paste(unique(hgnc_symbol_col[nzchar(hgnc_symbol_col)]), collapse = "; "),
    omim_gene = paste(unique(omim_gene_col[nzchar(omim_gene_col)]), collapse = "; "),
    hgnc_id = paste(unique(hgnc_id_col[nzchar(hgnc_id_col)]), collapse = "; "),
    ensembl_location = paste(unique(grch38_loc_col[nzchar(grch38_loc_col)]), collapse = "; "),
    ensembl_id = paste(unique(grch38_id_col[nzchar(grch38_id_col)]), collapse = "; ")
  ) |>
    mutate(
      search_text = paste(
        panel_name,
        disease_group,
        disease_sub_group,
        entity_names,
        phenotypes,
        mode_of_inheritance,
        mode_of_pathogenicity,
        publications,
        evidence,
        tags,
        alias,
        alias_name,
        gene_name,
        gene_symbol,
        hgnc_symbol,
        omim_gene,
        hgnc_id,
        ensembl_location,
        ensembl_id,
        sep = " ; "
      ),
      search_text_norm = normalise(search_text)
    )
}

build_ge_panel_index_from_ids <- function(panel_ids, api_key = NA_character_) {
  rows <- list()
  problems <- list()
  total_api_calls <- 0L
  
  for (pid in as.integer(panel_ids)) {
    got <- fetch_ge_panel_detail(pid, api_key = api_key)
    total_api_calls <- total_api_calls + as.integer(got$api_calls %||% 0L)
    
    if (!identical(got$status, "ok")) {
      problems[[length(problems) + 1L]] <- tibble(
        panel_id = pid,
        status = as.character(got$status),
        message = as.character(got$message %||% "")
      )
      next
    }
    
    row_i <- tryCatch(
      build_ge_single_panel_index(got$obj),
      error = function(e) {
        problems[[length(problems) + 1L]] <<- tibble(
          panel_id = pid,
          status = "no_genes",
          message = conditionMessage(e)
        )
        NULL
      }
    )
    
    if (!is.null(row_i)) {
      rows[[length(rows) + 1L]] <- row_i
    }
  }
  
  index_df <- if (length(rows) > 0) {
    bind_rows(rows) |>
      distinct(panel_id, .keep_all = TRUE)
  } else {
    empty_ge_index()
  }
  
  problems_df <- if (length(problems) > 0) {
    bind_rows(problems)
  } else {
    tibble(
      panel_id = integer(),
      status = character(),
      message = character()
    )
  }
  
  list(
    index = index_df,
    problems = problems_df,
    api_calls = total_api_calls
  )
}

run_panelapprex_end_to_end <- function(query, causal_gene, panel_subset_ids, path_core_rds) {
  t0 <- proc.time()[["elapsed"]]
  
  df_core <- readRDS(path_core_rds) |>
    as_tibble()
  
  df_core_subset <- df_core |>
    filter(as.integer(panel_id) %in% as.integer(panel_subset_ids))
  
  index_df <- collapse_panelapp_subset(df_core_subset) |>
    distinct(panel_id, .keep_all = TRUE)
  
  ranked <- rank_panels(index_df, query = query, causal_gene = causal_gene)$ranked
  
  elapsed_total_seconds <- proc.time()[["elapsed"]] - t0
  
  summarise_search_result(
    ranked_df = ranked,
    query = query,
    causal_gene = causal_gene,
    method = "PanelAppRex_local",
    subset_size = length(panel_subset_ids),
    elapsed_total_seconds = elapsed_total_seconds,
    n_api_calls = 0L,
    n_panels_requested = length(panel_subset_ids),
    n_panels_indexed = nrow(index_df)
  )
}

run_ge_end_to_end <- function(query, causal_gene, panel_subset_ids, api_key = NA_character_) {
  t0 <- proc.time()[["elapsed"]]
  
  ge_build <- build_ge_panel_index_from_ids(
    panel_ids = panel_subset_ids,
    api_key = api_key
  )
  
  ranked <- rank_panels(ge_build$index, query = query, causal_gene = causal_gene)$ranked
  
  elapsed_total_seconds <- proc.time()[["elapsed"]] - t0
  
  out <- summarise_search_result(
    ranked_df = ranked,
    query = query,
    causal_gene = causal_gene,
    method = "GE_PanelApp_API",
    subset_size = length(panel_subset_ids),
    elapsed_total_seconds = elapsed_total_seconds,
    n_api_calls = ge_build$api_calls,
    n_panels_requested = length(panel_subset_ids),
    n_panels_indexed = nrow(ge_build$index)
  )
  
  list(
    result = out,
    problems = ge_build$problems
  )
}

# ------------------------------------------------
# Inputs
# ------------------------------------------------

stop_if_missing(path_benchmark)
stop_if_missing(path_core_rds)

api_key <- read_api_key(path_api_key)

df_bench <- read.delim(
  path_benchmark,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

case_col <- extract_case_col(df_bench)
query_col <- extract_query_col(df_bench)
causal_col <- extract_causal_col(df_bench)

df_cases <- df_bench |>
  transmute(
    case_study = as.character(.data[[case_col]]),
    query = as.character(.data[[query_col]]),
    causal_gene = as.character(.data[[causal_col]])
  )

if (case_index < 1L || case_index > nrow(df_cases)) {
  stop("case_index is out of range for benchmark.tsv", call. = FALSE)
}

case_row <- df_cases[case_index, , drop = FALSE]

# ------------------------------------------------
# Build matched subset once
# ------------------------------------------------

meta_fetch <- fetch_all_ge_panel_metadata(api_key = api_key)

panel_subset_ids <- choose_panel_subset(
  meta_df = meta_fetch$df,
  seed_panels = seed_panels,
  n_random_extra_panels = n_random_extra_panels,
  random_seed = random_seed
)

results_tag_final <- sprintf(
  "cold_start_seed%s_plus%d_total%d_case%s",
  paste(seed_panels, collapse = "-"),
  n_random_extra_panels,
  length(panel_subset_ids),
  case_index
)

path_results_tsv <- file.path(path_data, sprintf("benchmark_results_%s.tsv", results_tag_final))
path_subset_tsv  <- file.path(path_data, sprintf("benchmark_subset_%s.tsv", results_tag_final))
path_summary_tsv <- file.path(path_data, sprintf("benchmark_summary_%s.tsv", results_tag_final))
path_problems_tsv <- file.path(path_data, sprintf("benchmark_ge_problems_%s.tsv", results_tag_final))

subset_df <- tibble(
  subset_size = length(panel_subset_ids),
  panel_id = as.integer(panel_subset_ids),
  seeded = as.integer(as.integer(panel_id) %in% as.integer(seed_panels))
)

write_tsv(subset_df, path_subset_tsv)

# ------------------------------------------------
# Run cold-start end-to-end benchmark
# ------------------------------------------------

panelapprex_res <- run_panelapprex_end_to_end(
  query = case_row$query[[1]],
  causal_gene = case_row$causal_gene[[1]],
  panel_subset_ids = panel_subset_ids,
  path_core_rds = path_core_rds
) |>
  mutate(case_study = case_row$case_study[[1]], .before = 1)

ge_run <- run_ge_end_to_end(
  query = case_row$query[[1]],
  causal_gene = case_row$causal_gene[[1]],
  panel_subset_ids = panel_subset_ids,
  api_key = api_key
)

ge_res <- ge_run$result |>
  mutate(case_study = case_row$case_study[[1]], .before = 1)

results_df <- bind_rows(panelapprex_res, ge_res) |>
  select(
    case_study,
    method,
    subset_size,
    query,
    query_tokens,
    n_query_tokens,
    elapsed_total_seconds,
    top_panel_id,
    top_panel_name,
    n_candidate_panels_returned,
    causal_gene,
    causal_gene_found_in_any_returned_panel,
    best_panel_contains_causal_gene,
    n_api_calls,
    n_panels_requested,
    n_panels_indexed
  ) |>
  arrange(method)

summary_df <- results_df |>
  mutate(
    seeded_panel_398_in_subset = as.integer(398L %in% panel_subset_ids),
    n_subset_panels = length(panel_subset_ids),
    subset_panel_ids = paste(panel_subset_ids, collapse = ", ")
  )

write_tsv(results_df, path_results_tsv)
write_tsv(summary_df, path_summary_tsv)
write_tsv(ge_run$problems, path_problems_tsv)

options(
  tibble.print_min = Inf,
  tibble.print_max = Inf,
  width = 1000
)

print(results_df, n = Inf, width = Inf)

if (nrow(ge_run$problems) > 0) {
  cat("\nGE problems:\n")
  print(ge_run$problems, n = Inf, width = Inf)
}

cat("\nSaved:\n")
cat("  results: ", path_results_tsv, "\n", sep = "")
cat("  subset:  ", path_subset_tsv, "\n", sep = "")
cat("  summary: ", path_summary_tsv, "\n", sep = "")
cat("  problems:", path_problems_tsv, "\n", sep = "")
cat("Subset size requested: ", length(panel_subset_ids), "\n", sep = "")
cat("Panel ids: ", paste(panel_subset_ids, collapse = ", "), "\n", sep = "")