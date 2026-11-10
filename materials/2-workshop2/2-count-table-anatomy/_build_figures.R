# Draw the microViz figures the Module 2 slides show at the end: the two
# composition bar charts and the two heatmaps.
#
# Run by hand when the data or this script changes, from the website root:
#   RENV_CONFIG_AUTOLOADER_ENABLED=FALSE \
#     Rscript materials/2-workshop2/2-count-table-anatomy/_build_figures.R
#
# Writes figures/*.png and results/top10_share.csv next to the slides. The
# slides include the pngs as images, so rendering the site never needs microViz
# or Bioconductor. The calls here are the ones printed on the slides, plus
# sizing and the community strip.

suppressPackageStartupMessages({library(tidyverse); library(microViz)})
here <- "materials/2-workshop2/2-count-table-anatomy"
out <- file.path(here, "figures"); dir.create(out, showWarnings = FALSE, recursive = TRUE)
source("materials/2-workshop2/R/count_tables.R")

d <- read_counts("instructional")
groups <- read_csv("materials/2-workshop2/4-differential/results/sample_groups.csv",
                   show_col_types = FALSE) |>
  select(sample_id, community)

species_ps <- as_counts(table = "microbiome", data = d) |> ps_join(groups, by = "sample_id")
genes_ps   <- as_counts(table = "rnaseq", data = d) |> ps_join(groups, by = "sample_id")

# --- composition bar charts ---------------------------------------------------
bar <- function(ps, rank, file) {
  p <- ps |>
    comp_barplot(tax_level = rank, n_taxa = 10, sample_order = "bray", label = NULL) +
    cowplot::theme_cowplot(16) +
    labs(y = "share of the sample") +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          axis.line.x = element_blank())
  ggsave(file.path(out, file), p, width = 10, height = 4.0, dpi = 160, bg = "white")
}
bar(species_ps, "Genus", "bar_species.png")
bar(genes_ps, "symbol", "bar_genes.png")

# --- heatmaps -----------------------------------------------------------------
community_col <- c(lactobacillus = "#2a9d8f", anaerobe = "#e76f51")
heat <- function(ps, rank, n, file) {
  ps_log <- ps |>
    tax_transform("compositional", rank = rank) |>
    tax_transform("log2", zero_replace = "halfmin", chain = TRUE)
  h <- ps_log |>
    comp_heatmap(
      taxa = tax_top(ps_log, n = n), name = "log2 share", grid_col = NA,
      sample_anno = sampleAnnotation(
        community = anno_sample_cat("community", col = community_col,
                                    legend_title = "community",
                                    size = grid::unit(7, "mm"))
      ),
      row_names_gp = grid::gpar(fontsize = if (n > 20) 9.5 else 11),
      heatmap_legend_param = list(labels_gp = grid::gpar(fontsize = 12),
                                  title_gp = grid::gpar(fontsize = 13, fontface = "bold"))
    )
  png(file.path(out, file), width = 10, height = 4.2, units = "in", res = 160)
  ComplexHeatmap::draw(h, merge_legends = TRUE)
  dev.off()
}
heat(species_ps, "Genus", 20, "heat_species.png")
heat(genes_ps, "symbol", 25, "heat_genes.png")

# --- how much of a sample the ten largest features hold -----------------------
top10 <- function(ps, rank, label) {
  ps |>
    tax_agg(rank) |>
    ps_melt() |>
    as_tibble() |>
    group_by(OTU) |>
    mutate(total = sum(Abundance)) |>
    ungroup() |>
    mutate(top10 = dense_rank(desc(total)) <= 10) |>
    group_by(Sample) |>
    summarise(share = sum(Abundance[top10]) / sum(Abundance)) |>
    summarise(table = label, features = n_distinct(phyloseq::taxa_names(tax_agg(ps, rank) |> ps_get())),
              median_share = median(share), min_share = min(share), max_share = max(share))
}
top10_share <- bind_rows(top10(species_ps, "Genus", "species"),
                         top10(genes_ps, "symbol", "genes"))
write_csv(top10_share, file.path(here, "results", "top10_share.csv"))
print(top10_share)
cat("module 2 figures written to", out, "\n")
