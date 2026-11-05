# Precompute the reveal figures for the Monday analog session.
#
# results/analog_cards.csv is copied from the private data generator, the
# printed card data joined to its bag-type key (data-gen/R/05_analog_cards.R:
# output/analog_cards/cards.csv + truth/analog_cards_truth.csv). The key is a
# spoiler, so the slides show it only after a warning slide.
#
# Run by hand when the cards change, from the website root:
#   RENV_CONFIG_AUTOLOADER_ENABLED=FALSE \
#     Rscript materials/2-workshop2/1-reads-to-counts/_build_results.R

suppressPackageStartupMessages({library(tidyverse); library(vegan); library(ggdendro)})
here <- "materials/2-workshop2/1-reads-to-counts"
out <- file.path(here, "results")

cards <- read_csv(file.path(out, "analog_cards.csv"), show_col_types = FALSE)
colours <- c("orange", "green", "blue", "pink", "red", "sky_blue", "yellow", "grey")
props <- as.matrix(cards[colours]) / cards$depth
rownames(props) <- cards$card_id

# Bray-Curtis on proportions, which is what the card fronts show.
bc <- vegdist(props, method = "bray")
pco <- cmdscale(bc, k = 3, eig = TRUE)
pct <- 100 * pco$eig[1:3] / sum(pco$eig[pco$eig > 0])
cards |>
  select(card_id, source, depth) |>
  mutate(PCo1 = pco$points[, 1], PCo2 = pco$points[, 2], PCo3 = pco$points[, 3]) |>
  write_csv(file.path(out, "analog_pcoa.csv"))
tibble(axis = paste0("PCo", 1:3), pct = round(pct, 1)) |>
  write_csv(file.path(out, "analog_pcoa_axes.csv"))

dd <- dendro_data(as.dendrogram(hclust(bc, method = "average")))
write_csv(dd$segments, file.path(out, "analog_dendro_segments.csv"))
dd$labels |> rename(card_id = label) |> write_csv(file.path(out, "analog_dendro_labels.csv"))
