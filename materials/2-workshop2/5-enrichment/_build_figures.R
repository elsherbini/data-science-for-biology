# Draw the microViz taxonomic-tree figure for the end of the Module 5 slides.
#
# Run by hand when the data or this script changes, from the website root:
#   RENV_CONFIG_AUTOLOADER_ENABLED=FALSE \
#     Rscript materials/2-workshop2/5-enrichment/_build_figures.R
#
# Writes figures/taxatree_outcome.png next to the slides. The slides include
# the png as an image, so rendering the site never needs microViz. The method
# follows https://david-barnett.github.io/microViz/articles/web-only/modelling-taxa.html

suppressPackageStartupMessages({library(tidyverse); library(microViz)})
here <- "materials/2-workshop2/5-enrichment"
out <- file.path(here, "figures"); dir.create(out, showWarnings = FALSE, recursive = TRUE)
source("materials/2-workshop2/R/count_tables.R")

d <- read_counts("instructional")

# --- outer-epithelium genes by community, for the "individual genes" slide ----
# Counts per million at week 18, so the slides can draw the boxplot without
# the full count table.
res_dir <- file.path(here, "results")
groups <- read_csv("materials/2-workshop2/4-differential/results/sample_groups.csv",
                   show_col_types = FALSE) |>
  filter(time_point == "week_18") |>
  select(sample_id, community)
cpm <- t(t(d$genes) / colSums(d$genes)) * 1e6
epi <- d$gene_info |> filter(symbol %in% c("FLG", "MUC22"))
cpm[epi$ensembl_gene_id, groups$sample_id, drop = FALSE] |>
  as_tibble(rownames = "ensembl_gene_id") |>
  pivot_longer(-ensembl_gene_id, names_to = "sample_id", values_to = "cpm") |>
  left_join(epi |> select(ensembl_gene_id, symbol), by = "ensembl_gene_id") |>
  left_join(groups, by = "sample_id") |>
  select(sample_id, community, symbol, cpm) |>
  write_csv(file.path(res_dir, "epithelial_genes_week18.csv"))

# Week 18 only, the same samples as Module 4's species test, with term as the
# reference so a positive estimate means higher in preterm.
ps <- as_counts(table = "microbiome", data = d) |>
  ps_filter(time_point == "week_18") |>
  ps_mutate(outcome = factor(outcome, c("term", "preterm")))

# One linear model per taxon at every rank, on log2 proportions. The
# prevalence filter keeps the 53 species seen in at least 10% of samples,
# which is the species set the slides rank.
models <- ps |>
  tax_prepend_ranks() |>
  tax_transform("compositional", rank = "species", keep_counts = TRUE) |>
  tax_filter(min_prevalence = 0.1, undetected = 0, use_counts = TRUE) |>
  taxatree_models(
    type = lm,
    trans = "log2", trans_args = list(zero_replace = "halfmin"),
    ranks = NULL,
    variables = "outcome",
    verbose = FALSE
  )
stats <- taxatree_models2stats(models)

p <- stats |>
  taxatree_label(rank == "Genus", p.value < 0.05 | prevalence > 0.7) |>
  taxatree_plots() |>
  pluck(1) |>
  taxatree_plot_labels(
    taxon_renamer = function(x) str_remove(x, "^[A-Za-z]: "),
    fun = ggrepel::geom_label_repel, x_nudge = 0.7, size = 3.6
  ) +
  labs(title = NULL)

ggsave(file.path(out, "taxatree_outcome.png"), p,
       width = 9, height = 6.5, dpi = 160, bg = "white")
