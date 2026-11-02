#  04_preprocess.R  (v2)  

SEED          <- 13
TEST_FRACTION <- 0.30
CLIP_Z        <- 3
USE_LOG1P     <- TRUE     
SCALE_TYPE    <- "iqr"    
COMPARE_OPTS  <- TRUE      
PCA_SAMPLE    <- 15000
PLOT_SAMPLE   <- 20000
SAVE_LEAKY    <- TRUE

OUT_DIR <- "data"; FIG_DIR <- "figures"; REP_DIR <- "reports"


for (p in c("data.table", "ggplot2", "scales"))
  if (!requireNamespace(p, quietly = TRUE))
    stop(sprintf("Install it first:  install.packages('%s')", p))

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(scales)
})
for (d in c(file.path(OUT_DIR, "processed"), FIG_DIR, REP_DIR))
  dir.create(d, showWarnings = FALSE, recursive = TRUE)

theme_ids <- function(base = 11) {
  theme_minimal(base_size = base) +
    theme(panel.grid.minor = element_blank(),
          plot.title    = element_text(face = "bold", size = base + 2),
          plot.subtitle = element_text(colour = "grey35", size = base - 1),
          plot.caption  = element_text(colour = "grey45", size = base - 2, hjust = 0),
          strip.text    = element_text(face = "bold"),
          legend.position = "bottom")
}
PAL_8 <- c(Benign = "#2E7D32", DDoS = "#C62828", DoS = "#EF6C00",
           Mirai = "#6A1B9A", Recon = "#1565C0", Spoofing = "#00838F",
           Web = "#AD1457", BruteForce = "#5D4037")
save_fig <- function(p, name, w = 9, h = 5.5) {
  path <- file.path(FIG_DIR, paste0(name, ".png"))
  ggsave(path, p, width = w, height = h, dpi = 200, bg = "white")
  cat("  figure ->", path, "\n"); invisible(path)
}


ws         <- readRDS(file.path(OUT_DIR, "interim", "working_set.rds"))
KEEP_FEATS <- readRDS(file.path(OUT_DIR, "interim", "keep_features.rds"))
cat(sprintf("Stage 4: %s rows x %d retained features\n",
            format(nrow(ws), big.mark = ","), length(KEEP_FEATS)))



set.seed(SEED)
idx_dt   <- data.table(row = seq_len(nrow(ws)), y = ws$label_8)
test_idx <- idx_dt[, .(row = sample(row, floor(TEST_FRACTION * .N))), by = y]$row
is_test  <- logical(nrow(ws)); is_test[test_idx] <- TRUE

cat(sprintf("train %s | test %s\n",
            format(sum(!is_test), big.mark = ","),
            format(sum(is_test),  big.mark = ",")))

strat <- ws[, .(n = .N), by = label_8]
strat <- merge(strat, ws[!is_test, .(n_train = .N), by = label_8], by = "label_8")
strat <- merge(strat, ws[ is_test, .(n_test  = .N), by = label_8], by = "label_8")
strat[, `:=`(p_train = round(100 * n_train / sum(n_train), 3),
             p_test  = round(100 * n_test  / sum(n_test),  3))]
strat[, max_dev := abs(p_train - p_test)]
setorder(strat, -n); print(strat)
cat("largest train/test proportion gap:", round(max(strat$max_dev), 4), "pp\n")
fwrite(strat, file.path(REP_DIR, "04_split_stratification.csv"))

y8_train    <- ws[!is_test, label_8];      y8_test    <- ws[is_test, label_8]
ybin_train  <- ws[!is_test, label_binary]; ybin_test  <- ws[is_test, label_binary]
yfine_train <- ws[!is_test, label_fine];   yfine_test <- ws[is_test, label_fine]

X_all <- as.matrix(ws[, ..KEEP_FEATS]); storage.mode(X_all) <- "double"
rm(ws, idx_dt); gc(verbose = FALSE)
cat(sprintf("feature matrix: %.2f GB\n", as.numeric(object.size(X_all)) / 1e9))



slog <- function(x) sign(x) * log1p(abs(x))

robust_spread <- function(x, m, type) {
  s <- switch(type,
              mad = mad(x, center = m, constant = 1.4826),
              iqr = unname(diff(quantile(x, c(0.25, 0.75), names = FALSE))) / 1.349,
              stop("SCALE_TYPE must be 'mad' or 'iqr'"))
  if (is.finite(s) && s > 0) return(list(s = s, src = toupper(type)))
  s <- sd(x)
  if (is.finite(s) && s > 0) return(list(s = s, src = sprintf("SD (%s = 0)", toupper(type))))
  list(s = 1, src = "1 (all spreads = 0)")
}

fit_transform <- function(X, rows = NULL) {
  if (is.null(rows)) rows <- seq_len(nrow(X))
  p <- ncol(X)
  raw_med <- centre <- scl <- numeric(p); src <- character(p)
  for (j in seq_len(p)) {
    x <- X[rows, j]; x <- x[is.finite(x)]
    if (!length(x)) { raw_med[j] <- 0; centre[j] <- 0; scl[j] <- 1
                      src[j] <- "empty"; next }
    raw_med[j] <- median(x)
    xt <- if (USE_LOG1P) slog(x) else x
    centre[j]  <- median(xt)
    rs <- robust_spread(xt, centre[j], SCALE_TYPE)
    scl[j] <- rs$s; src[j] <- rs$src
  }
  list(raw_median = raw_med, centre = centre, scale = scl,
       source = src, columns = colnames(X), log1p = USE_LOG1P)
}

apply_transform <- function(X, tf) {
  for (j in seq_len(ncol(X))) {
    x <- X[, j]
    x[!is.finite(x)] <- tf$raw_median[j]         
    if (tf$log1p) x <- slog(x)
    z <- (x - tf$centre[j]) / tf$scale[j]
    X[, j] <- pmax(pmin(z, CLIP_Z), -CLIP_Z)
  }
  X
}

cat(sprintf("\nfitting on TRAIN only  (log1p = %s, scale = %s)\n",
            USE_LOG1P, toupper(SCALE_TYPE)))
tr_rows <- which(!is_test)
tf_safe <- fit_transform(X_all, tr_rows); gc(verbose = FALSE)

tf_dt <- data.table(column = tf_safe$columns, raw_median = tf_safe$raw_median,
                    centre = tf_safe$centre, scale = tf_safe$scale,
                    source = tf_safe$source)
cat("\nspread statistic used:\n"); print(table(tf_dt$source))
print(tf_dt[order(source, column)], nrows = 40)
fwrite(tf_dt, file.path(REP_DIR, "04_transform_parameters.csv"))

## ---- clipping under the chosen transform ---------------------------------
clip_rate <- sapply(seq_along(KEEP_FEATS), function(j) {
  x <- X_all[tr_rows, j]; x <- x[is.finite(x)]
  if (USE_LOG1P) x <- slog(x)
  mean(abs((x - tf_safe$centre[j]) / tf_safe$scale[j]) > CLIP_Z)
})
clip_dt <- data.table(column = KEEP_FEATS, clip_share = clip_rate,
                      source = tf_safe$source)
setorder(clip_dt, -clip_share)
cat("\nclipping rate (train):\n"); print(clip_dt, nrows = 40)
cat(sprintf("mean %.2f%% | max %.2f%% (%s)\n",
            100 * mean(clip_rate), 100 * max(clip_rate),
            clip_dt[1, column]))
fwrite(clip_dt, file.path(REP_DIR, "04_clipping_rates.csv"))

## ---- apply, verify, save --------------------------------------------------
X_safe <- apply_transform(X_all, tf_safe)
stopifnot(all(is.finite(X_safe)))
stopifnot(max(abs(X_safe)) <= CLIP_Z + 1e-9)
cat("\npost-transform check: all finite, max |z| =", round(max(abs(X_safe)), 6), "\n")

saveRDS(list(X_train = X_safe[!is_test, , drop = FALSE],
             X_test  = X_safe[ is_test, , drop = FALSE],
             y8_train = y8_train,   y8_test = y8_test,
             ybin_train = ybin_train, ybin_test = ybin_test,
             yfine_train = yfine_train, yfine_test = yfine_test,
             features = KEEP_FEATS, transform = tf_safe,
             is_test = is_test, seed = SEED, test_fraction = TEST_FRACTION,
             clip_z = CLIP_Z),
        file.path(OUT_DIR, "processed", "split_safe.rds"))
cat("saved -> data/processed/split_safe.rds\n")



cat("\n*** LEAKY MODE: fitting on train+test. For comparison only. ***\n")
tf_leaky <- fit_transform(X_all, seq_len(nrow(X_all)))

cmp <- data.table(column = KEEP_FEATS,
                  centre_safe = tf_safe$centre, centre_leaky = tf_leaky$centre,
                  scale_safe  = tf_safe$scale,  scale_leaky  = tf_leaky$scale,
                  source_safe = tf_safe$source, source_leaky = tf_leaky$source)
cmp[, d_centre := abs(centre_leaky - centre_safe) / pmax(abs(centre_safe), 1e-12)]
cmp[, d_scale  := abs(scale_leaky  - scale_safe)  / pmax(abs(scale_safe),  1e-12)]
cmp[, source_changed := source_safe != source_leaky]
setorder(cmp, -d_scale)

cat("\nparameter divergence, safe vs leaky:\n")
print(cmp[, .(column, scale_safe = signif(scale_safe, 5),
              scale_leaky = signif(scale_leaky, 5),
              d_scale = signif(d_scale, 3), d_centre = signif(d_centre, 3),
              source_changed)], nrows = 40)
cat("\nstatistic changed in:", sum(cmp$source_changed), "columns\n")
cat("max relative scale divergence :", signif(max(cmp$d_scale), 4), "\n")
cat("max relative centre divergence:", signif(max(cmp$d_centre), 4), "\n")
fwrite(cmp, file.path(REP_DIR, "04_safe_vs_leaky_parameters.csv"))

if (SAVE_LEAKY) {
  X_leaky <- apply_transform(X_all, tf_leaky)
  saveRDS(list(X_train = X_leaky[!is_test, , drop = FALSE],
               X_test  = X_leaky[ is_test, , drop = FALSE],
               y8_train = y8_train, y8_test = y8_test,
               ybin_train = ybin_train, ybin_test = ybin_test,
               features = KEEP_FEATS, transform = tf_leaky,
               is_test = is_test, seed = SEED),
          file.path(OUT_DIR, "processed", "split_leaky.rds"))
  row_gap <- rowMeans(abs(X_safe - X_leaky))
  cat("saved -> data/processed/split_leaky.rds\n")
  cat("mean |safe - leaky| per row:", signif(mean(row_gap), 4),
      "| max:", signif(max(row_gap), 4), "\n")
  rm(X_leaky); gc(verbose = FALSE)
}

## Fig 12
p12 <- ggplot(tf_dt, aes(reorder(column, scale), scale, fill = source)) +
  geom_col() + coord_flip() +
  scale_y_log10(labels = comma) +
  labs(title = sprintf("Fig 12. Spread statistic used per feature (%s)", toupper(SCALE_TYPE)),
       subtitle = sprintf("%d of %d features fall back to SD because the robust spread is exactly zero -\na deliberate, reported behaviour for features whose modal value covers most rows.",
                          sum(tf_dt$source != toupper(SCALE_TYPE)), nrow(tf_dt)),
       x = NULL, y = "scale parameter (log10)", fill = NULL) +
  theme_ids()
save_fig(p12, "12_scale_statistic", h = 9)

## Fig 13
if (COMPARE_OPTS) {
  cat("\nrecomputing the scaling-option comparison for Fig 13 ...\n")
  clip_for <- function(logt, type, cz) {
    sapply(seq_len(ncol(X_all)), function(j) {
      x <- X_all[tr_rows, j]; x <- x[is.finite(x)]
      if (logt) x <- slog(x)
      m  <- median(x)
      rs <- robust_spread(x, m, type)
      mean(abs((x - m) / rs$s) > cz)
    })
  }
  opts <- list(`MAD`            = clip_for(FALSE, "mad", CLIP_Z),
               `IQR/1.349`      = clip_for(FALSE, "iqr", CLIP_Z),
               `log1p + MAD`    = clip_for(TRUE,  "mad", CLIP_Z),
               `log1p + IQR`    = clip_for(TRUE,  "iqr", CLIP_Z))
  odt <- rbindlist(lapply(names(opts), function(n)
    data.table(option = n, column = KEEP_FEATS, clip = opts[[n]])))
  odt[, option := factor(option, levels = names(opts))]
  fwrite(dcast(odt, column ~ option, value.var = "clip"),
         file.path(REP_DIR, "04_clipping_options.csv"))

  worst <- odt[option == "MAD"][order(-clip)][1:12, column]
  p13 <- ggplot(odt[column %in% worst],
                aes(reorder(column, clip), 100 * clip, fill = option)) +
    geom_col(position = "dodge") + coord_flip() +
    scale_fill_manual(values = c(`MAD` = "#C62828", `IQR/1.349` = "#EF6C00",
                                 `log1p + MAD` = "#6A1B9A", `log1p + IQR` = "#1565C0")) +
    labs(title = "Fig 13. Why the scaling method was changed",
         subtitle = "Share of training rows clipped at +/-3, twelve worst features under median/MAD.\nlog1p fixes heavy tails (Rate, Variance); IQR fixes bimodality (TCP, UDP). Neither alone suffices.",
         x = NULL, y = "% of training rows clipped", fill = NULL) +
    theme_ids()
  save_fig(p13, "13_clipping_comparison", h = 7)
  print(dcast(odt[column %in% worst], column ~ option, value.var = "clip")[
        , lapply(.SD, function(x) if (is.numeric(x)) round(100 * x, 2) else x)])
}

##Fig 14
set.seed(SEED)
pick   <- sample(tr_rows, min(PLOT_SAMPLE, length(tr_rows)))
n_uniq <- apply(X_all[pick, , drop = FALSE], 2, function(x) length(unique(x)))
show   <- KEEP_FEATS[order(-n_uniq)][1:min(8, length(KEEP_FEATS))]
jj     <- match(show, KEEP_FEATS)

before <- melt(as.data.table(X_all [pick, jj, drop = FALSE]),
               measure.vars = show, variable.name = "feature")[, stage := "Before"]
after  <- melt(as.data.table(X_safe[pick, jj, drop = FALSE]),
               measure.vars = show, variable.name = "feature")[, stage := "After"]
before[, value := slog(value)]
ba <- rbind(before, after)[is.finite(value)]
ba[, stage := factor(stage, levels = c("Before", "After"))]

p14 <- ggplot(ba, aes(value, fill = stage)) +
  geom_density(alpha = 0.6, colour = NA) +
  facet_wrap(~ feature + stage, scales = "free", ncol = 4) +
  scale_fill_manual(values = c(Before = "grey55", After = "#1565C0")) +
  labs(title = "Fig 14. Effect of log1p, robust scaling and clipping",
       subtitle = sprintf("Signed log1p, then median/%s z-score clipped to +/-%g, fitted on the training partition only.\n'Before' is on the signed log10 scale; 'After' is on the native z scale.",
                          toupper(SCALE_TYPE), CLIP_Z),
       x = NULL, y = "density", fill = NULL) +
  theme_ids(9)
save_fig(p14, "14_scaling_effect", w = 12, h = 9)

## Fig 15
set.seed(SEED)
pp  <- sample(tr_rows, min(PCA_SAMPLE, length(tr_rows)))
pca <- prcomp(X_safe[pp, , drop = FALSE], center = TRUE, scale. = FALSE)
ev  <- 100 * (pca$sdev^2) / sum(pca$sdev^2)
pdt <- data.table(PC1 = pca$x[, 1], PC2 = pca$x[, 2],
                  class = y8_train[match(pp, tr_rows)])
p15 <- ggplot(pdt, aes(PC1, PC2, colour = class)) +
  geom_point(alpha = 0.35, size = 0.7) +
  scale_colour_manual(values = PAL_8) +
  guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1))) +
  labs(title = "Fig 15. PCA of the preprocessed training data",
       subtitle = sprintf("PC1 %.1f%% / PC2 %.1f%% of variance. Diagnostic only - the paper selects features rather than projecting them.",
                          ev[1], ev[2]),
       x = sprintf("PC1 (%.1f%%)", ev[1]), y = sprintf("PC2 (%.1f%%)", ev[2]),
       colour = NULL) + theme_ids()
save_fig(p15, "15_pca", w = 9, h = 7)

scree <- data.table(PC = seq_along(ev), var = ev)[1:min(20, length(ev))]
scree[, cum := cumsum(var)]
p16 <- ggplot(scree, aes(PC)) +
  geom_col(aes(y = var), fill = "#1565C0") +
  geom_line(aes(y = cum), colour = "#C62828", linewidth = 0.8) +
  geom_point(aes(y = cum), colour = "#C62828") +
  geom_hline(yintercept = 90, linetype = "dashed", colour = "grey40") +
  labs(title = "Fig 16. Scree plot with cumulative variance",
       subtitle = sprintf("%d components reach 90%% of variance. Shown as the alternative we did NOT take.",
                          which(scree$cum >= 90)[1]),
       x = "Principal component", y = "% variance (bars) / cumulative % (line)") +
  theme_ids()
save_fig(p16, "16_scree")
fwrite(scree, file.path(REP_DIR, "04_pca_variance.csv"))

## Fig 17
p17 <- ggplot(cmp, aes(reorder(column, d_scale), 100 * d_scale,
                       fill = source_changed)) +
  geom_col() + coord_flip() +
  scale_fill_manual(values = c(`FALSE` = "grey60", `TRUE` = "#C62828"),
                    labels = c("same statistic", "statistic changed"), name = NULL) +
  labs(title = "Fig 17. Transform parameters: leakage-safe vs leaky",
       subtitle = sprintf("Relative change in the scale parameter when fitted on train+test instead of train alone.\nMaximum %.2f%% across %d features - small at this sample size, and reported as such.",
                          100 * max(cmp$d_scale), nrow(cmp)),
       x = NULL, y = "relative change in scale parameter (%)") +
  theme_ids()
save_fig(p17, "17_safe_vs_leaky", h = 9)


si <- c(
  sprintf("timestamp      : %s", format(Sys.time())),
  sprintf("R version      : %s", R.version.string),
  sprintf("platform       : %s", R.version$platform),
  sprintf("data.table     : %s", as.character(packageVersion("data.table"))),
  sprintf("seed           : %d", SEED),
  sprintf("TEST_FRACTION  : %.2f", TEST_FRACTION),
  sprintf("USE_LOG1P      : %s", USE_LOG1P),
  sprintf("SCALE_TYPE     : %s", SCALE_TYPE),
  sprintf("CLIP_Z         : %g", CLIP_Z),
  sprintf("features       : %d", length(KEEP_FEATS)),
  sprintf("train rows     : %d", sum(!is_test)),
  sprintf("test rows      : %d", sum(is_test)),
  sprintf("robust spread  : %d columns", sum(tf_dt$source == toupper(SCALE_TYPE))),
  sprintf("SD fallback    : %d columns", sum(tf_dt$source != toupper(SCALE_TYPE))),
  sprintf("mean clip rate : %.4f", mean(clip_rate)),
  sprintf("max clip rate  : %.4f (%s)", max(clip_rate), clip_dt[1, column]),
  sprintf("max scale divergence safe vs leaky: %.6g", max(cmp$d_scale))
)
writeLines(si, file.path(REP_DIR, "04_run_record.txt"))
cat("\n", paste(si, collapse = "\n"), "\n")
cat("\nStage 4 complete. -> data/processed/split_safe.rds\n")
