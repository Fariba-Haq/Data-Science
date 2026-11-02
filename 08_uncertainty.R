
SEED    <- 13
N_PERM  <- 2000       
ALPHA   <- 0.05
OUT_DIR <- "data"; FIG_DIR <- "figures"; REP_DIR <- "reports"

for (p in c("data.table", "ggplot2", "scales"))
  if (!requireNamespace(p, quietly = TRUE))
    stop(sprintf("Install it first:  install.packages('%s')", p))
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(scales)
})
for (d in c(FIG_DIR, REP_DIR)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

theme_ids <- function(base = 11) {
  theme_minimal(base_size = base) +
    theme(panel.grid.minor = element_blank(),
          plot.title    = element_text(face = "bold", size = base + 2),
          plot.subtitle = element_text(colour = "grey35", size = base - 1),
          plot.caption  = element_text(colour = "grey45", size = base - 2, hjust = 0),
          strip.text    = element_text(face = "bold"),
          legend.position = "bottom")
}
save_fig <- function(p, name, w = 9, h = 5.5) {
  path <- file.path(FIG_DIR, paste0(name, ".png"))
  ggsave(path, p, width = w, height = h, dpi = 200, bg = "white")
  cat("  figure ->", path, "\n"); invisible(path)
}

sel_obj <- readRDS(file.path(OUT_DIR, "processed", "selection.rds"))
diag    <- sel_obj$diag
FEATS   <- sel_obj$features
d       <- length(FEATS)
K_VALUES<- sel_obj$k_values
K_MAIN  <- sel_obj$k_main

ctrl_rows <- which(diag$config == "control_base")
configs   <- unique(diag$config)

cat(sprintf("Stage 8: %d draws, %d control draws, k = %s\n",
            nrow(diag), length(ctrl_rows), paste(K_VALUES, collapse = ", ")))


jaccard_matrix <- function(rank_mat, k) {
  M <- nrow(rank_mat)
  Z <- matrix(0L, M, d)
  for (i in seq_len(M)) Z[i, rank_mat[i, seq_len(k)]] <- 1L
  I <- tcrossprod(Z)                 # intersection sizes
  I / (2 * k - I)
}

JM <- list()
for (sel in c("L1", "RF")) {
  RM <- if (sel == "L1") sel_obj$rank_l1 else sel_obj$rank_rf
  for (k in K_VALUES) JM[[paste(sel, k)]] <- jaccard_matrix(RM, k)
}
cat("pairwise Jaccard matrices built\n")

cross_vals <- function(J, rows_a, rows_b, exclude_self = FALSE) {
  v <- J[rows_a, rows_b, drop = FALSE]
  if (exclude_self) v[upper.tri(v, diag = TRUE)] <- NA_real_
  as.numeric(v[!is.na(v)])
}

dist_dt <- list()
for (sel in c("L1", "RF")) for (k in K_VALUES) {
  J <- JM[[paste(sel, k)]]
  for (cf in configs) {
    rows <- which(diag$config == cf)
    is_ctrl <- cf == "control_base"
    v <- cross_vals(J, rows, ctrl_rows, exclude_self = is_ctrl)
    dist_dt[[length(dist_dt) + 1]] <- data.table(
      selector = sel, k = k, config = cf, arm = diag$arm[rows[1]],
      js = diag$js_realised_vs_base[rows[1]],
      n_pairs = length(v),
      mean = mean(v), sd = sd(v),
      lo = quantile(v, 0.025, names = FALSE),
      med = median(v),
      hi = quantile(v, 0.975, names = FALSE))
  }
}
dist_dt <- rbindlist(dist_dt)
fwrite(dist_dt, file.path(REP_DIR, "08_jaccard_distributions.csv"))

cat("\n--- NOISE FLOOR: control vs control ---\n")
print(dist_dt[config == "control_base",
              .(selector, k, n_pairs,
                mean = round(mean, 4), sd = round(sd, 4),
                ci = sprintf("[%.3f, %.3f]", lo, hi))])


set.seed(SEED)
perm_test <- function(J, rows_a, rows_ctrl, n_perm = N_PERM) {
  obs  <- mean(cross_vals(J, rows_a, rows_ctrl))
  pool <- c(rows_a, rows_ctrl)
  na   <- length(rows_a)
  null <- numeric(n_perm)
  for (b in seq_len(n_perm)) {
    p  <- sample(pool)
    null[b] <- mean(cross_vals(J, p[seq_len(na)], p[-seq_len(na)]))
  }
  list(obs = obs, null_mean = mean(null),
       p = (1 + sum(null <= obs)) / (n_perm + 1))
}

test_dt <- list()
for (sel in c("L1", "RF")) for (kk in K_VALUES) {
  J <- JM[[paste(sel, kk)]]
  floor_mean <- dist_dt[selector == sel & k == kk & config == "control_base", mean]
  stopifnot(length(floor_mean) == 1L)          # guard: must be a scalar
  for (cf in setdiff(configs, "control_base")) {
    rows <- which(diag$config == cf)
    tt <- perm_test(J, rows, ctrl_rows)
    test_dt[[length(test_dt) + 1]] <- data.table(
      selector = sel, k = kk, config = cf, arm = diag$arm[rows[1]],
      js = diag$js_realised_vs_base[rows[1]],
      obs_mean = tt$obs, noise_floor = floor_mean,
      drop_vs_floor = floor_mean - tt$obs, p_raw = tt$p)
  }
}
test_dt <- rbindlist(test_dt)

stopifnot(!anyDuplicated(test_dt[, .(selector, k, config)]))
test_dt[, p_holm := p.adjust(p_raw, method = "holm"), by = .(selector, k)]
test_dt[, significant := p_holm < ALPHA]
fwrite(test_dt, file.path(REP_DIR, "08_permutation_tests.csv"))

cat("\n", strrep("=", 78), "\n", sep = "")
cat(sprintf("  DROP-ONE ARM WITH INTERVALS  (k = %d)\n", K_MAIN))
cat(strrep("=", 78), "\n", sep = "")

do_dt <- merge(dist_dt[arm %in% c("drop_one", "control", "reference") & k == K_MAIN],
               test_dt[k == K_MAIN, .(selector, k, config, p_holm, significant)],
               by = c("selector", "k", "config"), all.x = TRUE)
setorder(do_dt, selector, -mean)
print(do_dt[, .(selector, config, n_pairs,
                mean = round(mean, 4),
                ci = sprintf("[%.3f, %.3f]", lo, hi),
                p_holm = ifelse(is.na(p_holm), NA, signif(p_holm, 3)),
                sig = ifelse(is.na(significant), "(floor)",
                             ifelse(significant, "YES", "no")))], nrows = 30)

cat("\n--- do the two selectors agree on the drop-one ORDERING? ---\n")
ord <- dcast(do_dt[arm == "drop_one"], config ~ selector, value.var = "mean")
if (all(c("L1", "RF") %in% names(ord))) {
  ord[, rank_L1 := rank(-L1)][, rank_RF := rank(-RF)]
  print(ord[order(rank_L1)])
  rho <- suppressWarnings(cor(ord$L1, ord$RF, method = "spearman"))
  cat(sprintf("\nSpearman rank correlation between selectors: %.3f\n", rho))

  cat(if (is.na(rho))
        "  -> Ordering could not be computed.\n"
      else if (rho <= -0.3)
        "  -> The orderings are ANTI-CORRELATED. The selectors disagree on which\n     category is most disruptive. Report the aggregate effect only; a\n     category ranking is not supported and must not be claimed.\n"
      else if (rho < 0.5)
        "  -> The per-category ordering does NOT replicate across selectors.\n     Report the aggregate effect; do NOT claim a category ranking.\n"
      else
        "  -> The ordering is consistent across selectors and can be reported.\n")
}


draw_dt <- list()
for (sel in c("L1", "RF")) for (k in K_VALUES) {
  J <- JM[[paste(sel, k)]]
  for (r in seq_len(nrow(diag))) {
    if (r %in% ctrl_rows) next
    draw_dt[[length(draw_dt) + 1]] <- data.table(
      selector = sel, k = k, arm = diag$arm[r], config = diag$config[r],
      js = diag$js_realised_vs_base[r],
      jaccard = mean(J[r, ctrl_rows]))
  }
}
draw_dt <- rbindlist(draw_dt)
fwrite(draw_dt, file.path(REP_DIR, "08_draw_level.csv"))

cat("\n--- dose-response: linear vs hyperbolic (k =", K_MAIN, ") ---\n")
fit_dt <- rbindlist(lapply(c("L1", "RF"), function(sel) {
  dd <- draw_dt[selector == sel & k == K_MAIN & js > 0]
  lin <- lm(jaccard ~ js, data = dd)
  hyp <- tryCatch(nls(jaccard ~ a / (1 + b * js), data = dd,
                      start = list(a = 1, b = 1)),
                  error = function(e) NULL)
  data.table(selector = sel, n = nrow(dd),
             slope_linear = round(coef(lin)[2], 4),
             r2_linear    = round(summary(lin)$r.squared, 4),
             aic_linear   = round(AIC(lin), 1),
             aic_hyperbolic = if (is.null(hyp)) NA_real_ else round(AIC(hyp), 1),
             better = if (is.null(hyp)) "linear (hyperbolic did not converge)"
                      else if (AIC(hyp) < AIC(lin) - 2) "HYPERBOLIC"
                      else if (AIC(lin) < AIC(hyp) - 2) "LINEAR"
                      else "indistinguishable")
}))
print(fit_dt)
fwrite(fit_dt, file.path(REP_DIR, "08_dose_response_fits.csv"))


## Fig 26: the drop-one arm, with intervals
p26 <- ggplot(do_dt, aes(reorder(config, mean), mean, colour = arm)) +
  geom_hline(data = do_dt[config == "control_base"],
             aes(yintercept = mean), linetype = "dashed", colour = "grey40") +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.25, linewidth = 0.5) +
  geom_point(size = 2.4) +
  coord_flip(ylim = c(0, 1)) +
  facet_wrap(~selector) +
  labs(title = sprintf("Fig 26. Drop-one arm with 95%% intervals (k = %d)", K_MAIN),
       subtitle = "Each point is the mean Jaccard over every (configuration draw x control draw) pair.\nDashed line is the control-vs-control noise floor. Overlapping intervals are NOT distinguishable.",
       x = NULL, y = "Jaccard against control draws", colour = NULL) +
  theme_ids()
save_fig(p26, "26_drop_one_intervals", w = 11, h = 6)

## Fig 27: distributions, not summaries
p27 <- ggplot(draw_dt[k == K_MAIN], aes(jaccard, fill = arm)) +
  geom_density(alpha = 0.5, colour = NA) +
  geom_vline(data = dist_dt[config == "control_base" & k == K_MAIN],
             aes(xintercept = mean), linetype = "dashed", colour = "grey30") +
  facet_wrap(~selector) +
  labs(title = sprintf("Fig 27. Distribution of subset overlap by arm (k = %d)", K_MAIN),
       subtitle = "Dashed line: the control-vs-control mean. Separation between the arms is the effect;\noverlap between them is what a single point estimate would have hidden.",
       x = "Jaccard against control draws", y = "density", fill = NULL) +
  theme_ids()
save_fig(p27, "27_overlap_distributions", w = 10, h = 5)

#Fig 28: dose-response with both fits
dd <- draw_dt[k == K_MAIN & js > 0]
p28 <- ggplot(dd, aes(js, jaccard)) +
  geom_point(aes(colour = arm), alpha = 0.45, size = 1.3) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE,
              colour = "#1565C0", linewidth = 0.7, linetype = "dashed") +
  geom_smooth(method = "nls", formula = y ~ a / (1 + b * x), se = FALSE,
              method.args = list(start = list(a = 1, b = 1)),
              colour = "#C62828", linewidth = 0.9) +
  facet_wrap(~selector) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(title = sprintf("Fig 28. Dose-response, fitted at draw level (k = %d)", K_MAIN),
       subtitle = "Blue dashed: linear. Red solid: hyperbolic, the pre-registered prediction.\nSee reports/08_dose_response_fits.csv for the AIC comparison.",
       x = "Jensen-Shannon divergence from the base mix",
       y = "Jaccard against control draws", colour = NULL) +
  theme_ids()
save_fig(p28, "28_dose_response_fitted", w = 10, h = 5)

## Fig 29: does the effect hold at every k
p29 <- ggplot(dist_dt[arm != "control"], aes(factor(k), mean, colour = arm)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.15,
                                             dodge.width = 0.6),
             alpha = 0.55, size = 1.3) +
  geom_point(data = dist_dt[config == "control_base"],
             aes(factor(k), mean), colour = "black", size = 3, shape = 18) +
  facet_wrap(~selector) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(title = "Fig 29. Does the composition effect survive the choice of k?",
       subtitle = "Black diamonds: the control noise floor at each k. If the arms sit below it at every k,\nthe finding is not an artefact of the feature budget.",
       x = "top-k budget", y = "mean Jaccard against control draws", colour = NULL) +
  theme_ids()
save_fig(p29, "29_effect_across_k", w = 10, h = 5)


cat("\n", strrep("=", 78), "\n", sep = "")
cat("  WHAT CAN BE CLAIMED\n")
cat(strrep("=", 78), "\n", sep = "")

agg <- dist_dt[k == K_MAIN, .(configs = .N, mean_jaccard = round(mean(mean), 4),
                              lo = round(min(lo), 3), hi = round(max(hi), 3)),
               by = .(selector, arm)]
print(agg)

nsig <- test_dt[k == K_MAIN & arm == "drop_one",
                .(significant = sum(significant), of = .N), by = selector]
cat("\ndrop-one configurations significantly below the noise floor (Holm, k =",
    K_MAIN, "):\n"); print(nsig)

saveRDS(list(distributions = dist_dt, tests = test_dt, draws = draw_dt,
             fits = fit_dt, n_perm = N_PERM, seed = SEED),
        file.path(OUT_DIR, "processed", "uncertainty.rds"))

si <- c(
  sprintf("timestamp      : %s", format(Sys.time())),
  sprintf("R version      : %s", R.version.string),
  sprintf("seed           : %d", SEED),
  sprintf("permutations   : %d", N_PERM),
  sprintf("noise floor L1 (k=%d): %.4f", K_MAIN,
          dist_dt[selector == "L1" & k == K_MAIN & config == "control_base", mean]),
  sprintf("noise floor RF (k=%d): %.4f", K_MAIN,
          dist_dt[selector == "RF" & k == K_MAIN & config == "control_base", mean])
)
writeLines(si, file.path(REP_DIR, "08_run_record.txt"))
cat("\n", paste(si, collapse = "\n"), "\n")
cat("\nStage 8 complete. -> data/processed/uncertainty.rds\n")
