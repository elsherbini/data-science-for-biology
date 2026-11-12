# Precompute the Module 6 integration results.
#
# Run by hand from the website root, after modules 4 and 5:
#   RENV_CONFIG_AUTOLOADER_ENABLED=FALSE \
#     Rscript materials/2-workshop2/6-integration/_build_results.R
#
# One number per person from each table, then the three models the Friday
# discussion walks through. Nothing exotic: this is the simple baseline every
# fancier integration method has to beat.

suppressPackageStartupMessages({library(tidyverse); library(msigdbr); library(pROC)})
here <- "materials/2-workshop2/6-integration"
out <- file.path(here, "results"); dir.create(out, showWarnings = FALSE, recursive = TRUE)
source("materials/2-workshop2/R/count_tables.R")

d <- read_counts("instructional")
groups <- read_csv("materials/2-workshop2/4-differential/results/sample_groups.csv",
                   show_col_types = FALSE)

# Host score: mean z-scored log-CPM of two inflammatory Hallmark sets.
hallmark <- msigdbr(species = "Homo sapiens", collection = "H")
inflam <- hallmark |>
  filter(gs_name %in% c("HALLMARK_INFLAMMATORY_RESPONSE", "HALLMARK_TNFA_SIGNALING_VIA_NFKB")) |>
  pull(ensembl_gene) |> unique()

cpm <- sweep(d$genes, 2, colSums(d$genes), "/") * 1e6
expressed <- rownames(cpm)[rowMeans(cpm) >= 5]
lcpm <- log2(cpm[expressed, ] + 1)
set_genes <- intersect(inflam, rownames(lcpm))
cat("inflammatory genes used for the score:", length(set_genes), "of", length(inflam), "\n")
z <- t(scale(t(lcpm[set_genes, ])))
host_score <- colMeans(z, na.rm = TRUE)

per_sample <- groups |>
  select(sample_id, participant_id, time_point, outcome, community, lacto_fraction) |>
  mutate(anaerobe_fraction = 1 - lacto_fraction,
         host_score = host_score[sample_id])
write_csv(per_sample, file.path(out, "scores_per_sample.csv"))

per_person <- per_sample |>
  group_by(participant_id, outcome) |>
  summarise(micro_score = mean(anaerobe_fraction),
            host_score = mean(host_score), .groups = "drop") |>
  mutate(preterm = as.integer(outcome == "preterm"),
         micro_z = as.numeric(scale(micro_score)),
         host_z = as.numeric(scale(host_score)))
write_csv(per_person, file.path(out, "scores_per_person.csv"))

fits <- list(
  `microbiome only` = glm(preterm ~ micro_z, per_person, family = binomial),
  `host only`       = glm(preterm ~ host_z, per_person, family = binomial),
  `both`            = glm(preterm ~ micro_z + host_z, per_person, family = binomial)
)

model_tbl <- imap_dfr(fits, function(f, nm) {
  r <- roc(per_person$preterm, predict(f, type = "link"), quiet = TRUE)
  ci <- as.numeric(ci.auc(r))
  tibble(model = nm, auc = as.numeric(auc(r)), auc_lo = ci[1], auc_hi = ci[3],
         aic = AIC(f))
})
lrt <- anova(fits$`microbiome only`, fits$both, test = "Chisq")
model_tbl$lrt_p_vs_micro <- c(NA, NA, lrt$`Pr(>Chi)`[2])
write_csv(model_tbl, file.path(out, "model_comparison.csv"))
print(model_tbl)

coef_tbl <- imap_dfr(fits, function(f, nm) {
  broom_like <- summary(f)$coefficients
  as_tibble(broom_like, rownames = "term") |>
    mutate(model = nm, odds_ratio = exp(Estimate)) |>
    rename(estimate = Estimate, std_error = `Std. Error`, p_value = `Pr(>|z|)`)
})
write_csv(coef_tbl, file.path(out, "model_coefficients.csv"))

roc_tbl <- imap_dfr(fits, function(f, nm) {
  r <- roc(per_person$preterm, predict(f, type = "link"), quiet = TRUE)
  tibble(model = nm, fpr = 1 - r$specificities, tpr = r$sensitivities)
})
write_csv(roc_tbl, file.path(out, "roc_curves.csv"))

cat("correlation between the two scores:",
    sprintf("%.3f", cor(per_person$micro_z, per_person$host_z)), "\n")
cat("module 6 results written to", out, "\n")

# ---------------------------------------------------------------------------
# Predictive is not causal. Rank every per-person predictor by AUC and label
# it with its role in the DAG that generated the outcome (data-gen: preterm
# depends on mean anaerobe fraction and a latent host-inflammation flag, and
# on nothing else).
# ---------------------------------------------------------------------------

cyto <- read_csv(file.path(dataset_path("instructional"), "06_cytokines.csv"),
                 show_col_types = FALSE)
per_person_all <- d$metadata |>
  left_join(cyto, by = "sample_id") |>
  group_by(participant_id) |>
  summarise(nugent = mean(nugent, na.rm = TRUE), ph = mean(ph, na.rm = TRUE),
            clue_cells = mean(clue_cells, na.rm = TRUE),
            il_1b = mean(log10(il_1b + 0.01)), tnfa = mean(log10(tnfa + 0.01)),
            il_8 = mean(log10(il_8 + 0.01)),
            batch_b2 = mean(rnaseq_batch == "B2"), age = first(age)) |>
  inner_join(per_person, by = "participant_id")

predictors <- tribble(
  ~var,          ~label,                         ~role,
  "micro_score", "anaerobe fraction (species)",  "cause",
  "host_score",  "inflammatory score (genes)",   "cause",
  "nugent",      "Nugent score",                 "readout of a cause",
  "ph",          "vaginal pH",                   "readout of a cause",
  "clue_cells",  "clue cells",                   "readout of a cause",
  "il_1b",       "IL-1β (cytokine)",        "bystander",
  "tnfa",        "TNF-α (cytokine)",        "bystander",
  "il_8",        "IL-8 (cytokine)",              "bystander",
  "batch_b2",    "RNA sequencing batch",         "technical",
  "age",         "age",                          "unrelated"
)
predictor_ranking <- predictors |>
  mutate(auc = map_dbl(var, ~ as.numeric(auc(roc(per_person_all$preterm, per_person_all[[.x]],
                                                 quiet = TRUE, direction = "auto")))),
         cor_with_micro = map_dbl(var, ~ cor(per_person_all[[.x]], per_person_all$micro_score,
                                            use = "complete.obs")))
write_csv(predictor_ranking, file.path(out, "predictor_ranking.csv"))
print(predictor_ranking)

# ---------------------------------------------------------------------------
# mixOmics on one row per person: each feature averaged over the three
# visits. That avoids repeated measures, and it matches how the outcome was
# generated (from the mean anaerobe fraction across visits; week 18 alone
# gives the species table an AUC of only 0.57).
# Species: CLR, present in at least 10% of samples. Genes: log-CPM, the 2,000
# most variable expressed genes. Both z-scored by mixOmics (scale = TRUE).
# ---------------------------------------------------------------------------

suppressPackageStartupMessages(library(mixOmics))
select <- dplyr::select   # mixOmics masks both
map <- purrr::map
set.seed(20261113)

people <- per_person |>
  left_join(per_sample |> filter(time_point == "week_18") |>
              select(participant_id, sample_id), by = "participant_id") |>
  left_join(groups |> filter(time_point == "week_18") |>
              select(participant_id, community), by = "participant_id")
by_person <- function(m) {   # features x samples -> people x features, mean over visits
  pid <- d$metadata$participant_id[match(colnames(m), d$metadata$sample_id)]
  t(sapply(people$participant_id, function(p) rowMeans(m[, pid == p, drop = FALSE])))
}
sp <- d$species[rowMeans(d$species > 0) >= 0.10, ]
clr <- function(m) { l <- log(m + 0.5); t(t(l) - colMeans(l)) }
X_sp <- by_person(clr(sp))
X_gn_all <- by_person(lcpm)
X_gn <- X_gn_all[, order(-apply(X_gn_all, 2, var))[1:2000]]
Y <- factor(people$outcome, levels = c("term", "preterm"))
md18 <- people |> mutate(lacto_fraction = 1 - micro_score)   # names used below
sym <- deframe(d$gene_info |> select(ensembl_gene_id, symbol))
cat("person-level blocks: species", ncol(X_sp), " genes", ncol(X_gn), " people", nrow(X_sp), "\n")

# --- sPLS, species x genes, no outcome -------------------------------------
keep_sp <- c(10, 10); keep_gn <- c(50, 50)
spls_fit <- spls(X_sp, X_gn, ncomp = 2, keepX = keep_sp, keepY = keep_gn, mode = "canonical")
spls_scores <- tibble(participant_id = people$participant_id, outcome = people$outcome,
                      community = people$community,
                      species_comp1 = spls_fit$variates$X[, 1],
                      genes_comp1 = spls_fit$variates$Y[, 1])
write_csv(spls_scores, file.path(out, "spls_scores.csv"))

spls_loadings <- bind_rows(
  tibble(table = "species", feature = colnames(X_sp), loading = spls_fit$loadings$X[, 1]),
  tibble(table = "genes", feature = colnames(X_gn), loading = spls_fit$loadings$Y[, 1])
) |>
  filter(loading != 0) |>
  mutate(label = if_else(table == "genes", coalesce(sym[feature], feature), feature))
write_csv(spls_loadings, file.path(out, "spls_loadings.csv"))

# How tightly the selected genes move together, against random genes.
sel_genes <- spls_loadings |> filter(table == "genes") |> pull(feature)
mean_offdiag <- function(m) { r <- cor(m); mean(abs(r[upper.tri(r)])) }
module_cor <- tibble(
  set = c("50 genes sPLS picked", "50 random genes"),
  mean_abs_cor = c(mean_offdiag(X_gn[, sel_genes]),
                   mean(replicate(50, mean_offdiag(X_gn[, sample(ncol(X_gn), 50)]))))
)
write_csv(module_cor, file.path(out, "spls_module_cor.csv"))
print(module_cor)

# Rerun on bootstrap resamples: how often is each gene picked on component 1?
n_boot <- 50
boot_sel <- map(seq_len(n_boot), function(b) {
  i <- sample(nrow(X_sp), replace = TRUE)
  xs <- X_sp[i, ]; xg <- X_gn[i, ]
  rownames(xs) <- rownames(xg) <- make.unique(rownames(xs))   # mixOmics wants unique rows
  f <- spls(xs, xg, ncomp = 1, keepX = keep_sp[1], keepY = keep_gn[1],
            mode = "canonical")
  colnames(X_gn)[f$loadings$Y[, 1] != 0]
})
spls_stability <- tibble(feature = unlist(boot_sel)) |>
  count(feature, name = "times_picked") |>
  mutate(freq = times_picked / n_boot,
         label = coalesce(sym[feature], feature),
         in_full_fit = feature %in% sel_genes) |>
  arrange(desc(freq))
write_csv(spls_stability, file.path(out, "spls_stability.csv"))
cat("genes ever picked in", n_boot, "bootstraps:", nrow(spls_stability),
    "; picked >= 90% of the time:", sum(spls_stability$freq >= 0.9), "\n")

# --- cross-validation helpers ----------------------------------------------
# Stratified folds, repeated. Every model sees exactly the same folds.
make_folds <- function(y, k = 5, repeats = 10, seed = 1) {
  set.seed(seed)
  map(seq_len(repeats), function(r) {
    f <- integer(length(y))
    for (lv in levels(y)) { i <- which(y == lv); f[i] <- sample(rep_len(seq_len(k), length(i))) }
    f
  })
}
folds <- make_folds(Y, k = 5, repeats = 10)
auc_of_scores <- function(y, s) as.numeric(auc(roc(y, s, quiet = TRUE, levels = c("term", "preterm"), direction = "<")))

cv_diablo <- function(X, y, design, keepX, ncomp = 2, folds) {
  map_dfr(seq_along(folds), function(r) {
    f <- folds[[r]]
    pred <- map(seq_len(max(f)), function(k) {
      tr <- f != k; te <- f == k
      fit <- block.splsda(map(X, ~ .x[tr, , drop = FALSE]), y[tr], ncomp = ncomp,
                          keepX = keepX, design = design)
      p <- predict(fit, newdata = map(X, ~ .x[te, , drop = FALSE]))
      tibble(i = which(te),
             combined = p$AveragedPredict[, "preterm", ncomp],
             weighted = p$WeightedPredict[, "preterm", ncomp],
             !!!map(set_names(names(X)), ~ p$predict[[.x]][, "preterm", ncomp]))
    }) |> bind_rows() |> arrange(i)
    cols <- setdiff(names(pred), "i")
    tibble(repeat_id = r, model = cols,
           auc = map_dbl(cols, ~ auc_of_scores(y[pred$i], pred[[.x]])))
  })
}
design_mat <- function(c) matrix(c(0, c, c, 0), 2, dimnames = list(c("species", "genes"), c("species", "genes")))
blocks <- list(species = X_sp, genes = X_gn)
keepX <- list(species = c(10, 10), genes = c(50, 50))

# --- the two-score baseline, cross-validated on the same folds -------------
base_df <- tibble(y = Y, micro = people$micro_score, host = people$host_score)
cv_glm <- map_dfr(seq_along(folds), function(r) {
  f <- folds[[r]]
  s <- numeric(nrow(base_df))
  for (k in seq_len(max(f))) {
    fit <- glm(y ~ micro + host, base_df[f != k, ], family = binomial)
    s[f == k] <- predict(fit, base_df[f == k, ])
  }
  tibble(repeat_id = r, model = "two scores, one glm", auc = auc_of_scores(Y, s))
})

# --- DIABLO at design 0.1 ---------------------------------------------------
cv_d <- cv_diablo(blocks, Y, design_mat(0.1), keepX, folds = folds)
diablo_cv <- bind_rows(cv_glm, cv_d) |>
  mutate(model = recode(model, species = "DIABLO, species table's vote",
                        genes = "DIABLO, gene table's vote",
                        combined = "DIABLO, both tables averaged",
                        weighted = "DIABLO, both tables weighted")) |>
  group_by(model) |>
  summarise(auc_sd = sd(auc), auc = mean(auc), .groups = "drop")
write_csv(diablo_cv, file.path(out, "diablo_cv.csv"))
print(diablo_cv)

# Full-data DIABLO fit for loadings and the cross-table correlation.
diablo_fit <- block.splsda(blocks, Y, ncomp = 2, keepX = keepX, design = design_mat(0.1))
diablo_loadings <- bind_rows(
  tibble(table = "species", feature = colnames(X_sp), loading = diablo_fit$loadings$species[, 1]),
  tibble(table = "genes", feature = colnames(X_gn), loading = diablo_fit$loadings$genes[, 1])
) |>
  filter(loading != 0) |>
  mutate(label = if_else(table == "genes", coalesce(sym[feature], feature), feature))
# Which group has the higher mean, the way plotLoadings colours its bars.
higher_in <- function(m, f) if_else(mean(m[Y == "preterm", f]) > mean(m[Y == "term", f]), "preterm", "term")
diablo_loadings <- diablo_loadings |>
  mutate(higher_in = map2_chr(table, feature, ~ higher_in(blocks[[.x]], .y)))
write_csv(diablo_loadings, file.path(out, "diablo_loadings.csv"))

# --- the design setting as a trade-off --------------------------------------
diablo_design <- map_dfr(c(0, 0.1, 0.5, 0.9), function(cc) {
  fit <- block.splsda(blocks, Y, ncomp = 2, keepX = keepX, design = design_mat(cc))
  cv <- cv_diablo(blocks, Y, design_mat(cc), keepX, folds = folds[1:5])
  tibble(design = cc,
         auc = mean(cv$auc[cv$model == "combined"]),
         cross_table_cor = abs(cor(fit$variates$species[, 1], fit$variates$genes[, 1])))
})
write_csv(diablo_design, file.path(out, "diablo_design.csv"))
print(diablo_design)

# --- shuffled labels: the training plot separates anyway --------------------
Y_shuf <- sample(Y)
shuffle_scores <- map_dfr(list(real = Y, shuffled = Y_shuf), function(yy) {
  f <- splsda(X_gn, yy, ncomp = 2, keepX = c(50, 50))
  tibble(comp1 = f$variates$X[, 1], comp2 = f$variates$X[, 2], label = yy)
}, .id = "labels")
write_csv(shuffle_scores, file.path(out, "shuffle_scores.csv"))

cv_splsda <- function(yy, folds) map_dbl(seq_along(folds), function(r) {
  f <- folds[[r]]; s <- numeric(length(yy))
  for (k in seq_len(max(f))) {
    fit <- splsda(X_gn[f != k, ], yy[f != k], ncomp = 2, keepX = c(50, 50))
    s[f == k] <- predict(fit, X_gn[f == k, , drop = FALSE])$predict[, "preterm", 2]
  }
  auc_of_scores(yy, s)
})
shuffle_cv <- bind_rows(
  tibble(labels = "real", auc = cv_splsda(Y, folds[1:5])),
  map_dfr(1:5, function(s) {
    set.seed(100 + s); ys <- sample(Y)
    tibble(labels = "shuffled", auc = cv_splsda(ys, make_folds(ys, repeats = 2, seed = s)))
  })
) |>
  group_by(labels) |>
  summarise(auc_sd = sd(auc), auc = mean(auc), .groups = "drop")
write_csv(shuffle_cv, file.path(out, "shuffle_cv.csv"))
print(shuffle_cv)

cat("integration results written to", out, "\n")
