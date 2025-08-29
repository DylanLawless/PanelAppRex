# In this script we use df_core to get the coordindates for each gene in panel_id, then make one bed file per panel.

# Rds format core data ----
path_data <- "../data"
path_PanelAppData_genes_combined_Rds <- paste0(path_data, "/PanelAppData_genes_combined_Rds")
df_core <- readRDS(file= path_PanelAppData_genes_combined_Rds)

# make_gene_locations_grch38.R
suppressPackageStartupMessages({
  library(GenomicRanges)
  library(GenomeInfoDb)
  library(ensembldb)
  # pick one EnsDb matching GRCh38. Use v86 if you want older Ensembl 86 gene bounds.
  library(EnsDb.Hsapiens.v86)
})

# Return a data.frame with one row per unique HGNC symbol:
# columns: symbol, chrom, start, end, location ("chr:start-end")
get_gene_locations_from_symbols <- function(symbols, flank = 0L, use_v86 = TRUE) {
  stopifnot(is.character(symbols))
  syms <- unique(toupper(trimws(symbols)))
  syms <- syms[!is.na(syms) & nzchar(syms)]
  if (!length(syms)) {
    return(data.frame(symbol = character(), chrom = character(), start = integer(), end = integer(), location = character()))
  }
  
  edb <- if (use_v86) EnsDb.Hsapiens.v86 else EnsDb.Hsapiens.v105
  
  gr <- genes(edb, filter = SymbolFilter(syms))
  
  if (length(gr) == 0L) {
    return(data.frame(symbol = syms, chrom = NA_character_, start = NA_integer_, end = NA_integer_, location = NA_character_))
  }
  
  # standardise chromosomes
  suppressWarnings(seqlevelsStyle(gr) <- "UCSC")
  std_chr <- paste0("chr", c(1:22, "X", "Y", "M"))
  gr <- gr[as.character(seqnames(gr)) %in% std_chr]
  
  if (length(gr) == 0L) {
    return(data.frame(symbol = syms, chrom = NA_character_, start = NA_integer_, end = NA_integer_, location = NA_character_))
  }
  
  mcols(gr)$gene_name <- toupper(as.character(mcols(gr)$gene_name))
  mcols(gr)$gene_biotype <- as.character(mcols(gr)$gene_biotype)
  
  # collapse to one interval per symbol on a single chromosome:
  # prefer protein_coding if multiple loci, else pick chromosome with widest covered width
  split_by_sym <- split(gr, mcols(gr)$gene_name)
  
  pick_one <- function(g) {
    if (length(g) == 1L) return(g)
    g_pc <- g[mcols(g)$gene_biotype == "protein_coding"]
    g_use <- if (length(g_pc)) g_pc else g
    chr <- as.character(seqnames(g_use))
    total_by_chr <- tapply(width(g_use), chr, sum)
    best_chr <- names(total_by_chr)[which.max(total_by_chr)]
    g_chr <- g_use[as.character(seqnames(g_use)) == best_chr]
    GRanges(seqnames = best_chr,
            ranges = IRanges(start = min(start(g_chr)), end = max(end(g_chr))),
            strand = "*",
            gene_name = mcols(g_chr)$gene_name[1])
  }
  
  picked <- suppressWarnings(unlist(GRangesList(lapply(split_by_sym, pick_one)), use.names = FALSE))
  if (length(picked) == 0L) {
    return(data.frame(symbol = syms, chrom = NA_character_, start = NA_integer_, end = NA_integer_, location = NA_character_))
  }
  
  # add flank
  flank <- as.integer(flank)
  if (flank > 0L) {
    start(picked) <- pmax(start(picked) - flank, 1L)
    end(picked) <- end(picked) + flank
  }
  
  df <- data.frame(
    symbol = as.character(mcols(picked)$gene_name),
    chrom  = as.character(seqnames(picked)),
    start  = as.integer(start(picked)),
    end    = as.integer(end(picked)),
    stringsAsFactors = FALSE
  )
  df$location <- paste0(df$chrom, ":", df$start, "-", df$end)
  
  # ensure one row per requested symbol and in input order
  df <- df[!duplicated(df$symbol), ]
  missing_syms <- setdiff(syms, df$symbol)
  if (length(missing_syms)) {
    df <- rbind(df, data.frame(symbol = missing_syms, chrom = NA_character_, start = NA_integer_, end = NA_integer_, location = NA_character_))
  }
  df[match(syms, df$symbol), , drop = FALSE]
}

# Get positions ----
syms <- unique(na.omit(df_core$gene_data.hgnc_symbol))
# syms <- head(syms, 200) # subset test
loc_map <- get_gene_locations_from_symbols(syms, flank = 0L, use_v86 = TRUE)
idx <- match(toupper(df_core$gene_data.hgnc_symbol), loc_map$symbol)
df_core$gene_location_grch38 <- loc_map$location[idx]

# validate manually ----
df_core <- df_core |> dplyr::select(gene_data.ensembl_genes.GRch38.90.location, gene_location_grch38, everything())
df_core |> dplyr::filter(gene_data.hgnc_symbol == "ACTA2") |> dplyr::select(gene_data.ensembl_genes.GRch38.90.location, gene_location_grch38, gene_data.hgnc_symbol)

# save one BED per panel_id from df_core$gene_location_grch38 ----
# ensure output directory exists
outdir <- file.path(path_data, "bed")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# safe filename helper
safe_name <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x <- gsub("[^A-Za-z0-9._-]+", "_", x)
  x[nchar(x) == 0L] <- "unknown_panel"
  x
}

# parse locations already present in df_core
loc <- df_core$gene_location_grch38
keep <- !is.na(loc) & nzchar(loc) & grepl(":", loc)

chrom <- sub(":.*", "", loc[keep])
start1 <- suppressWarnings(as.integer(sub(".*:(\\d+)-.*", "\\1", loc[keep])))
end1   <- suppressWarnings(as.integer(sub(".*-(\\d+)$", "\\1", loc[keep])))

df_parsed <- data.frame(
  panel_id = df_core$panel_id[keep],
  symbol = toupper(as.character(df_core$gene_data.hgnc_symbol[keep])),
  chrom = chrom,
  chromStart = pmax(start1 - 1L, 0L),
  chromEnd = end1,
  stringsAsFactors = FALSE
)

# keep only standard chromosomes and valid rows
std_chr <- paste0("chr", c(1:22, "X", "Y", "M"))
df_parsed <- df_parsed[!is.na(df_parsed$panel_id) & !is.na(df_parsed$chrom) & df_parsed$chrom %in% std_chr, ]
df_parsed$chrom <- factor(df_parsed$chrom, levels = std_chr)

# write one BED per panel_id
panels <- unique(df_parsed$panel_id)
invisible(lapply(panels, function(pid) {
  d <- df_parsed[df_parsed$panel_id == pid, c("chrom", "chromStart", "chromEnd", "symbol")]
  if (!nrow(d)) return(NULL)
  # ensure one row per gene symbol
  d <- d[!is.na(d$symbol) & nzchar(d$symbol), ]
  d <- d[!duplicated(d$symbol), ]
  if (!nrow(d)) return(NULL)
  # sort
  o <- order(d$chrom, d$chromStart, d$chromEnd, d$symbol)
  d <- d[o, ]
  names(d) <- c("chrom", "chromStart", "chromEnd", "name")
  outfile <- file.path(outdir, paste0("panel_", safe_name(pid), ".bed"))
  write.table(d, outfile, sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
  NULL
}))
