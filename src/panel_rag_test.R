# panel_rag_test.R

library(dplyr)
library(readr)
library(purrr)
library(stringr)
library(httr)
library(jsonlite)

# -----------------------------
# Config
# -----------------------------

path_data <- "../data"
path_panels_rds <- file.path(path_data, "PanelAppData_genes_combined_Rds")

# UniProt discussion table from your previous pipeline
# Columns assumed: SYMBOL, Protein.names,
#   Gene.Ontology..molecular.function., Function..CC.
path_uniprot_discussion <- file.path(
  path_data,
  "get_discussion",
  "df_report_discussion_.tsv"
)

# Output for RAG panel summaries
out_dir <- file.path(path_data, "rag_panels_test")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Test panel id (adjust as needed)
test_panel_id <- 398L

# OpenAI key and model
api_key <- readLines(
  "/Users/dylanlawless/Desktop/switzerlandomics_api_key.txt",
  warn = FALSE
) |> str_trim()

openai_model <- "gpt-4.1-mini"

# -----------------------------
# Load data
# -----------------------------

df_panels <- readRDS(path_panels_rds)

# Use entity_name as SYMBOL
panel_genes <- df_panels |>
  filter(panel_id == test_panel_id) |>
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

df_discussion <- read_tsv(path_uniprot_discussion, show_col_types = FALSE)

panel_discussion <- panel_genes |>
  left_join(df_discussion, by = "SYMBOL")

# -----------------------------
# Build compact gene context
# -----------------------------

make_gene_snippet <- function(sym, prot, go_mf, func_cc) {
  fields <- c(
    paste0("Gene: ", sym),
    if (!is.na(prot) && nzchar(prot)) paste0("Protein: ", prot) else NA_character_,
    if (!is.na(go_mf) && nzchar(go_mf)) paste0("Molecular function: ", go_mf) else NA_character_,
    if (!is.na(func_cc) && nzchar(func_cc)) paste0("Cellular function: ", func_cc) else NA_character_
  )
  fields <- fields[!is.na(fields)]
  paste(fields, collapse = " | ")
}

panel_discussion <- panel_discussion |>
  mutate(
    gene_snippet = pmap_chr(
      list(
        SYMBOL,
        Protein.names,
        `Gene.Ontology..molecular.function.`,
        `Function..CC.`
      ),
      make_gene_snippet
    )
  )

# Optionally limit number of genes for the test
max_genes <- 50L
panel_discussion_small <- panel_discussion |>
  slice_head(n = max_genes)

gene_context_text <- panel_discussion_small$gene_snippet |>
  paste(collapse = "\n")

# -----------------------------
# Build prompt for OpenAI
# -----------------------------

panel_meta <- panel_discussion_small |>
  distinct(
    panel_id,
    panel_name,
    disease_group,
    disease_sub_group,
    mode_of_inheritance
  ) |>
  slice(1)

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

spec_text <- '
You are enriching a single gene panel for better natural language search.
Use the metadata and gene context below to infer what this panel is really about.

Return a single JSON object with this exact structure:

{
  "panel_id": <integer>,
  "panel_name": "<string>",
  "summary": "<one or two sentences, specific to this panel>",
  "semantic_tags": ["<string>", "..."],
  "mechanism_tags": ["<string>", "..."],
  "gene_level_tags": ["<string>", "..."],
  "example_queries": ["<string>", "..."]
}

Rules:
- Only include terms that are highly characteristic of this panel.
- Prefer disease names, phenotypes, mechanisms, and gene specific concepts.
- Avoid generic words like "genetic", "disease", "syndrome", "panel", "mutation".
- Limit:
  - semantic_tags to at most 20 items
  - mechanism_tags to at most 10 items
  - gene_level_tags to at most 10 items
  - example_queries to at most 5 items
- Use short phrases, not full sentences, for tags.
- summary must be concise and suitable as a description on a search page.
'

full_user_content <- paste0(
  spec_text,
  "\n\n",
  "Panel metadata:\n\n",
  meta_text,
  "\n\nGene context (one line per gene):\n\n",
  gene_context_text
)

body <- list(
  model = openai_model,
  response_format = list(type = "json_object"),
  messages = list(
    list(
      role = "system",
      content = "You are a clinical genomics expert who writes precise, technically correct summaries."
    ),
    list(
      role = "user",
      content = full_user_content
    )
  )
)

# -----------------------------
# Call OpenAI API
# -----------------------------

res <- POST(
  url = "https://api.openai.com/v1/chat/completions",
  add_headers(
    `Content-Type` = "application/json",
    Authorization = paste("Bearer", api_key)
  ),
  body = toJSON(body, auto_unbox = TRUE)
)

stop_for_status(res)

txt <- content(res, as = "text", encoding = "UTF-8")
parsed <- fromJSON(txt, simplifyVector = TRUE)

rag_json <- parsed$choices[[1]]$message$content

rag_obj <- fromJSON(rag_json, simplifyVector = TRUE)

# -----------------------------
# Save outputs
# -----------------------------

rds_path <- file.path(
  out_dir,
  paste0("panel_", test_panel_id, "_rag.rds")
)
json_path <- file.path(
  out_dir,
  paste0("panel_", test_panel_id, "_rag.json")
)

saveRDS(rag_obj, rds_path)
write_json(rag_obj, json_path, auto_unbox = TRUE, pretty = TRUE)

cat("Saved RAG summary for panel", test_panel_id, "to:\n  ", rds_path, "\n  ", json_path, "\n")
