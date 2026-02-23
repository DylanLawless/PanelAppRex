# build_uniprotkb_up000005640_whole_proteome.R

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(stringr)
  library(readr)
  library(tidyr)
  library(jsonlite)
  library(tibble)
})

# ------------------------------------------------
# Config
# ------------------------------------------------

panelapprex_root  <- "~/mnt/atlas_data_big/data/panelapprex"

uniprot_dir <- file.path(panelapprex_root, "uniprot_up000005640")
uniprot_dat <- file.path(uniprot_dir, "UP000005640_9606.dat.gz")

out_dir_processed <- "~/mnt/atlas_data_big/data/panelapprex/uniprot_kb_up000005640_processed"
# out_dir_processed <- "../data/uniprot_kb_up000005640"
dir.create(out_dir_processed, showWarnings = FALSE, recursive = TRUE)

# Accession-level outputs
out_accession_rds <- file.path(out_dir_processed, "uniprot_accession_level.rds")
out_accession_tsv <- file.path(out_dir_processed, "uniprot_accession_level.tsv")

# CC topics long table (easy to explore, grep, QC)
out_cc_long_tsv <- file.path(out_dir_processed, "uniprot_cc_topics_long.tsv")

# Gene-level outputs
out_gene_rds <- file.path(out_dir_processed, "uniprot_gene_level.rds")
out_gene_tsv <- file.path(out_dir_processed, "uniprot_gene_level.tsv")

# Tier 1 topics (embedding/search layer)
tier1_topics <- c(
  "FUNCTION",
  "DISEASE"
  # "SUBCELLULAR LOCATION",
  # "SIMILARITY"
)

# ------------------------------------------------
# Helpers
# ------------------------------------------------

first_or_na <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(NA_character_)
  x[[1]]
}

trim_or_na <- function(x) {
  if (is.na(x)) return(NA_character_)
  x <- str_squish(x)
  if (!nzchar(x)) return(NA_character_)
  x
}

extract_accession <- function(ac_lines) {
  if (length(ac_lines) == 0) return(NA_character_)
  ac_string <- ac_lines |>
    str_remove("^AC\\s+") |>
    str_c(collapse = " ")
  parts <- str_split(ac_string, ";", simplify = TRUE) |>
    as.character() |>
    str_trim()
  parts <- parts[parts != ""]
  if (length(parts) == 0) return(NA_character_)
  parts[[1]]
}

extract_symbol <- function(gn_lines) {
  if (length(gn_lines) == 0) return(NA_character_)
  sym <- gn_lines |>
    str_c(collapse = " ") |>
    str_extract("Name=([^;]+)") |>
    str_remove("^Name=")
  trim_or_na(sym)
}

extract_protein_name <- function(de_lines) {
  if (length(de_lines) == 0) return(NA_character_)
  recname <- de_lines |> str_subset("RecName: Full=")
  if (length(recname) == 0) return(NA_character_)
  nm <- recname |>
    str_c(collapse = " ") |>
    str_extract("RecName: Full=([^;]+);") |>
    str_remove("^RecName: Full=") |>
    str_remove(";$")
  trim_or_na(nm)
}

parse_cc_topics <- function(cc_lines) {
  cc_df <- tibble(topic = character(), content = character())
  if (length(cc_lines) == 0) return(cc_df)
  
  cc_text <- cc_lines |>
    str_remove("^CC\\s+") |>
    str_c(collapse = " ") |>
    str_squish()
  
  if (!nzchar(cc_text)) return(cc_df)
  
  blocks <- str_extract_all(
    cc_text,
    "-!-\\s+[A-Z][A-Z ]+:\\s+.*?(?=\\s+-!-\\s+|$)"
  )[[1]]
  
  if (length(blocks) == 0) return(cc_df)
  
  map_dfr(blocks, function(block) {
    topic <- str_extract(block, "-!-\\s+[A-Z][A-Z ]+:") |>
      str_remove("^-!-\\s+") |>
      str_remove(":$") |>
      str_squish()
    
    content <- block |>
      str_remove("^-!-\\s+[A-Z][A-Z ]+:\\s+") |>
      str_squish()
    
    tibble(topic = topic, content = content)
  }) |>
    filter(!is.na(topic), nzchar(topic), !is.na(content), nzchar(content))
}

build_tier1_text <- function(cc_df) {
  if (nrow(cc_df) == 0) return(NA_character_)
  txt <- cc_df |>
    filter(topic %in% tier1_topics) |>
    transmute(piece = paste0(topic, ": ", content)) |>
    pull(piece) |>
    unique() |>
    paste(collapse = " || ") |>
    str_squish()
  if (!nzchar(txt)) return(NA_character_)
  txt
}

# Streaming record reader: avoids holding all lines + records in memory at once
parse_uniprot_dat_stream <- function(path, progress_every = 1000L) {
  con <- gzfile(path, open = "rt")
  on.exit(close(con), add = TRUE)
  
  results <- vector("list", 0)
  cc_long <- vector("list", 0)
  
  current <- character()
  rec_i <- 0L
  
  repeat {
    chunk <- readLines(con, n = 50000L, warn = FALSE)
    if (length(chunk) == 0) break
    
    for (ln in chunk) {
      current <- c(current, ln)
      
      if (identical(ln, "//")) {
        rec_i <- rec_i + 1L
        
        ac_lines <- current[str_starts(current, "AC  ")]
        de_lines <- current[str_starts(current, "DE  ")]
        gn_lines <- current[str_starts(current, "GN  ")]
        cc_lines <- current[str_starts(current, "CC  ")]
        
        accession <- extract_accession(ac_lines)
        symbol <- extract_symbol(gn_lines)
        protein_name <- extract_protein_name(de_lines)
        
        cc_df <- parse_cc_topics(cc_lines)
        tier1_text <- build_tier1_text(cc_df)
        
        results[[rec_i]] <- tibble(
          accession = accession,
          SYMBOL = symbol,
          Protein_name = protein_name,
          Tier1_text = tier1_text,
          cc_topics = list(cc_df)
        )
        
        if (nrow(cc_df) > 0 && !is.na(accession)) {
          cc_long[[rec_i]] <- cc_df |>
            mutate(accession = accession, .before = 1L)
        } else {
          cc_long[[rec_i]] <- tibble(accession = character(), topic = character(), content = character())
        }
        
        if (rec_i %% progress_every == 0L) {
          cat("Parsed records:", rec_i, "\n")
        }
        
        current <- character()
      }
    }
  }
  
  df_accession <- bind_rows(results)
  df_cc_long <- bind_rows(cc_long) |>
    filter(!is.na(accession), nzchar(accession))
  
  list(accession = df_accession, cc_long = df_cc_long)
}

# Extract arrays from cc_topics by topic, preserving all blocks
topic_array <- function(cc_df, topic_name) {
  if (is.null(cc_df) || nrow(cc_df) == 0) return(character())
  cc_df |>
    filter(topic == topic_name) |>
    pull(content) |>
    unique() |>
    str_squish() |>
    (\(x) x[nzchar(x)])()
}

# ------------------------------------------------
# Parse whole proteome (accession-level)
# ------------------------------------------------

stopifnot(file.exists(uniprot_dat))

cat("Parsing UniProt dat (streaming):", uniprot_dat, "\n")

parsed <- parse_uniprot_dat_stream(uniprot_dat, progress_every = 1000L)
df_uniprot_accession <- parsed$accession
df_cc_long <- parsed$cc_long

cat("Accession rows:", nrow(df_uniprot_accession), "\n")
cat("Distinct accessions:", n_distinct(df_uniprot_accession$accession), "\n")
cat("Distinct symbols:", n_distinct(df_uniprot_accession$SYMBOL), "\n")
cat("CC long rows:", nrow(df_cc_long), "\n")

# Save accession-level TSV (flat, portable)
write_tsv(
  df_uniprot_accession |>
    select(accession, SYMBOL, Protein_name, Tier1_text),
  out_accession_tsv
)

# Save accession-level RDS (keeps list-column cc_topics)
saveRDS(df_uniprot_accession, out_accession_rds)

# Save CC long table (portable exploration)
write_tsv(df_cc_long, out_cc_long_tsv)

cat("Saved:\n")
cat("  ", out_accession_tsv, "\n")
cat("  ", out_accession_rds, "\n")
cat("  ", out_cc_long_tsv, "\n")

cat("Aggregating to gene-level\n")

df_gene_level <- df_uniprot_accession |>
  dplyr::filter(!is.na(SYMBOL), nzchar(SYMBOL)) |>
  dplyr::group_by(SYMBOL) |>
  dplyr::summarise(
    protein_names = list(sort(unique(stats::na.omit(Protein_name)))),
    function_terms = list(sort(unique(unlist(purrr::map(cc_topics, topic_array, "FUNCTION"))))),
    disease_terms = list(sort(unique(unlist(purrr::map(cc_topics, topic_array, "DISEASE"))))),
    subcellular_location_terms = list(sort(unique(unlist(purrr::map(cc_topics, topic_array, "SUBCELLULAR LOCATION"))))),
    similarity_terms = list(sort(unique(unlist(purrr::map(cc_topics, topic_array, "SIMILARITY"))))),
    tier1_text = paste(sort(unique(stats::na.omit(Tier1_text))), collapse = " || "),
    .groups = "drop"
  ) |>
  dplyr::mutate(
    tier1_text = dplyr::na_if(stringr::str_squish(tier1_text), ""),
    json = purrr::pmap_chr(
      list(
        SYMBOL,
        protein_names,
        function_terms,
        disease_terms,
        subcellular_location_terms,
        similarity_terms,
        tier1_text
      ),
      function(symbol,
               protein_names,
               function_terms,
               disease_terms,
               subcellular_location_terms,
               similarity_terms,
               tier1_text) {
        
        obj <- list(
          symbol = symbol,
          protein_names = protein_names,
          "function" = function_terms,
          disease = disease_terms,
          subcellular_location = subcellular_location_terms,
          similarity = similarity_terms,
          tier1_text = tier1_text
        )
        
        jsonlite::toJSON(
          obj,
          auto_unbox = TRUE,
          pretty = FALSE,
          null = "null"
        )
      }
    )
  )

cat("Gene-level rows:", nrow(df_gene_level), "\n")

saveRDS(df_gene_level, out_gene_rds)

readr::write_tsv(
  df_gene_level |>
    dplyr::transmute(
      SYMBOL,
      tier1_text,
      json
    ),
  out_gene_tsv
)

cat("Saved:\n")
cat("  ", out_gene_rds, "\n")
cat("  ", out_gene_tsv, "\n")

cat("\n--- RAG1 test (gene-level) ---\n")

rag1 <- df_gene_level |>
  dplyr::filter(SYMBOL == "RAG1")

if (nrow(rag1) == 0) {
  cat("RAG1 not found in gene-level table.\n")
} else {
  cat("RAG1 tier1_text:\n")
  cat(rag1$tier1_text[[1]], "\n\n")
  
  cat("RAG1 JSON:\n")
  cat(rag1$json[[1]], "\n\n")
  
  cat("Counts:\n")
  cat("  protein_names:", length(rag1$protein_names[[1]]), "\n")
  cat("  function_terms:", length(rag1$function_terms[[1]]), "\n")
  cat("  disease_terms:", length(rag1$disease_terms[[1]]), "\n")
  cat("  subcellular_location_terms:", length(rag1$subcellular_location_terms[[1]]), "\n")
  cat("  similarity_terms:", length(rag1$similarity_terms[[1]]), "\n\n")
  
  cat("First few disease items:\n")
  cat(" - ", utils::head(rag1$disease_terms[[1]], 10), sep = "\n - ")
  cat("\n")
}

cat("\n--- end RAG1 test ---\n")

# ------------------------------------------------
# Note: PanelApp subsetting is intentionally not done here.
# Later we will join df_gene_level to PanelApp genes/panels/pathways and filter then.
# ------------------------------------------------
