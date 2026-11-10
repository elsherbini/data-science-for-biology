# count_tables.R ---------------------------------------------------------------
#
# One job: hand you a workshop count table in whichever shape you need next.
#
#   read_counts()                     plain matrices and tibbles
#   as_counts(table = "microbiome")   a phyloseq object microViz can pipe
#   as_counts(table = "rnaseq")       the same, built from the gene table
#
# The second one is the point of the week. microViz was written for microbiome
# data, so it expects a "taxonomy" table: one row per feature, columns running
# from coarse labels to fine ones. A gene annotation is the same thing wearing
# different words, so we hand microViz `locus_group` and `symbol` as the two
# ranks and every microViz verb then works on genes with no further coaxing.
#
#   library(microViz)
#   source("R/count_tables.R")
#
#   as_counts(table = "rnaseq") |>
#     tax_transform("clr", rank = "symbol") |>
#     dist_calc("euclidean") |>
#     ord_calc() |>
#     ord_plot(colour = "outcome")

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tibble)
})

# Folder names on disk, keyed by the short name you type.
WORKSHOP_DATASETS <- c(
  instructional = "instructional_dataset",
  group_a       = "group_a_probiotic",
  group_b       = "group_b_contraception"
)

# Where the data might be. The Codespace puts it in data/2-workshop2/, the
# website in datasets/2-workshop2/. Set WORKSHOP_DATA to override both.
#
# Note the "2-workshop2" in both. Workshop 1 also has a folder called
# data/instructional_dataset, holding completely different files, so a bare
# "data" root would find the wrong one.
dataset_path <- function(dataset = "instructional") {
  folder <- if (dataset %in% names(WORKSHOP_DATASETS)) WORKSHOP_DATASETS[[dataset]] else dataset
  roots <- c(Sys.getenv("WORKSHOP_DATA", unset = NA), "data/2-workshop2", "datasets/2-workshop2")
  if (requireNamespace("here", quietly = TRUE)) {
    roots <- c(roots, here::here("data/2-workshop2"), here::here("datasets/2-workshop2"))
  }
  for (r in roots[!is.na(roots)]) {
    p <- file.path(r, folder)
    if (dir.exists(p)) return(p)
  }
  stop("Could not find dataset '", dataset, "'. Looked for a folder called '",
       folder, "' under: ", paste(roots[!is.na(roots)], collapse = ", "),
       ". Set WORKSHOP_DATA to the folder that holds it.", call. = FALSE)
}

# A features-in-rows integer matrix from a csv whose first column is the ID.
read_count_matrix <- function(file) {
  d <- read_csv(file, show_col_types = FALSE, progress = FALSE)
  m <- as.matrix(d[, -1])
  rownames(m) <- d[[1]]
  storage.mode(m) <- "integer"
  m
}

#' Read one dataset.
#'
#' Returns a list: metadata, species, taxonomy, genes, gene_info, cytokines.
#' Both count matrices are features in rows, samples in columns, with columns
#' in the same order as metadata$sample_id.
read_counts <- function(dataset = "instructional") {
  p <- dataset_path(dataset)
  metadata <- read_csv(file.path(p, "01_metadata.csv"), show_col_types = FALSE, progress = FALSE)
  species  <- read_count_matrix(file.path(p, "02_microbiome_counts.csv"))
  genes    <- read_count_matrix(file.path(p, "04_rnaseq_counts.csv.gz"))
  stopifnot(setequal(metadata$sample_id, colnames(species)),
            setequal(metadata$sample_id, colnames(genes)))
  list(
    metadata  = metadata,
    species   = species[, metadata$sample_id, drop = FALSE],
    taxonomy  = read_csv(file.path(p, "03_microbiome_taxonomy.csv"), show_col_types = FALSE, progress = FALSE),
    genes     = genes[, metadata$sample_id, drop = FALSE],
    gene_info = read_csv(file.path(p, "05_gene_annotation.csv"), show_col_types = FALSE, progress = FALSE),
    cytokines = read_csv(file.path(p, "06_cytokines.csv"), show_col_types = FALSE, progress = FALSE),
    path      = p
  )
}

# Feature table -> character matrix of ranks, coarse to fine, rownames = feature id.
rank_matrix <- function(info, id_col, ranks) {
  ranks <- intersect(ranks, names(info))
  m <- as.matrix(info[, ranks, drop = FALSE])
  storage.mode(m) <- "character"
  rownames(m) <- info[[id_col]]
  m
}

#' Build a phyloseq object from either count table.
#'
#' @param dataset "instructional", "group_a", or "group_b".
#' @param table   "microbiome" for the species table, "rnaseq" for the gene table.
#' @param min_count,min_samples drop features seen fewer than `min_count` times
#'   in fewer than `min_samples` samples. Defaults keep everything.
#' @param fix     run microViz::tax_fix() to fill blank rank labels. Leave TRUE.
as_counts <- function(dataset = "instructional",
                      table = c("microbiome", "rnaseq"),
                      min_count = 0, min_samples = 0, fix = TRUE,
                      data = NULL) {
  table <- match.arg(table)
  if (!requireNamespace("phyloseq", quietly = TRUE)) {
    stop("as_counts() needs the phyloseq package.", call. = FALSE)
  }
  d <- data %||% read_counts(dataset)

  if (table == "microbiome") {
    counts <- d$species
    ranks  <- rank_matrix(d$taxonomy, "taxon",
                          c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "species"))
  } else {
    counts <- d$genes
    # A few genes have no symbol. Give them their Ensembl id, so every gene
    # keeps a label of its own.
    info <- d$gene_info |>
      mutate(symbol      = coalesce(symbol, ensembl_gene_id),
             locus_group = coalesce(locus_group, "unannotated"))
    ranks  <- rank_matrix(info, "ensembl_gene_id", c("locus_group", "symbol"))
  }

  # Features present in the annotation but not the counts, and vice versa.
  keep <- intersect(rownames(counts), rownames(ranks))
  counts <- counts[keep, , drop = FALSE]
  ranks  <- ranks[keep, , drop = FALSE]

  if (min_count > 0 || min_samples > 0) {
    enough <- rowSums(counts >= min_count) >= min_samples
    counts <- counts[enough, , drop = FALSE]
    ranks  <- ranks[enough, , drop = FALSE]
  }

  meta <- as.data.frame(d$metadata)
  rownames(meta) <- meta$sample_id

  ps <- phyloseq::phyloseq(
    phyloseq::otu_table(counts, taxa_are_rows = TRUE),
    phyloseq::sample_data(meta),
    phyloseq::tax_table(ranks)
  )
  if (fix && requireNamespace("microViz", quietly = TRUE)) {
    # tax_fix() treats labels shorter than 4 characters as unknown, which suits
    # taxonomy and is wrong for genes: B2M, FTL and about 470 other short
    # symbols would be pooled into one feature. Keep every gene symbol.
    ps <- microViz::tax_fix(ps, min_length = if (table == "rnaseq") 1 else 4,
                            verbose = FALSE)
  }
  ps
}

`%||%` <- function(x, y) if (is.null(x)) y else x
