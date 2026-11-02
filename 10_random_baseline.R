
#  10_random_baseline.R 



## CONFIG


SEED        <- 13
N_DRAWS     <- 10     
N_RANDOM    <- 20      
NUM_TREES   <- 100
TEST_SAMPLE <- 40000

OUT_DIR <- "data"; FIG_DIR <- "figures"; REP_DIR <- "reports"


need <- c("data.table", "ggplot2", "scales", "ranger")
miss <- need[!sapply(need, requireNamespace, quietly = TRUE)]
if (length(miss))
  stop("Install first: install.packages(c(",
       paste(sprintf("'%s'", miss), collapse = ", "), "))")
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(scales); library(ranger)
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

sp      <- readRDS(file.path(OUT_DIR, "processed", "split_safe.rds"))
sel_obj <- readRDS(file.path(OUT_DIR, "processed", "selection.rds"))
arms    <- readRDS(file.path(OUT_DIR, "processed", "arms.rds"))

Xtr <- sp$X_train; y8tr <- sp$y8_train
Xte <- sp$X_test;  y8te <- sp$y8_test
FEATS <- sp$features; d <- length(FEATS)
K <- sel_obj$k_main
safe_names <- make.names(FEATS, unique = TRUE)

diag  <- sel_obj$diag
draws <- arms$draws[diag$draw_key]
ctrl_rows <- which(diag$config == "control_base")[seq_len(N_DRAWS)]

set.seed(SEED)
te_idx <- unlist(lapply(unique(y8te), function(cl) {
  w <- which(y8te == cl)
  sample(w, max(1, round(length(w) * TEST_SAMPLE / length(y8te))))
}))
Xte_s <- Xte[te_idx, , drop = FALSE]; yte_s <- factor(y8te[te_idx])
LEVS <- levels(yte_s)

n_fits <- N_DRAWS * (1 + 2 + N_RANDOM + 2)
cat(sprintf("Stage 10: %d control draws x (1 all + 2 selected + %d random + 2 worst)\n",
            N_DRAWS, N_RANDOM))
cat(sprintf("          = %d model fits | test subsample %s rows | k = %d\n",
            n_fits, format(length(te_idx), big.mark = ","), K))

macro_f1 <- function(truth, pred, levs) {
  truth <- factor(truth, levels = levs); pred <- factor(pred, levels = levs)
  mean(sapply(levs, function(cl) {
    tp <- sum(pred == cl & truth == cl)
    fp <- sum(pred == cl & truth != cl)
    fn <- sum(pred != cl & truth == cl)
    if (tp == 0) return(0)
    p <- tp/(tp+fp); r <- tp/(tp+fn)
    if (p + r == 0) 0 else 2*p*r/(p+r)
  }))
}

fit_eval <- function(train_rows, cols, seed) {
  x <- Xtr[train_rows, cols, drop = FALSE]; colnames(x) <- safe_names[cols]
  rf <- ranger(x = x, y = factor(y8tr[train_rows], levels = LEVS),
               num.trees = NUM_TREES, num.threads = 0,
               seed = seed, verbose = FALSE)
  xt <- Xte_s[, cols, drop = FALSE]; colnames(xt) <- safe_names[cols]
  macro_f1(yte_s, predict(rf, data = xt, num.threads = 0)$predictions, LEVS)
}

set.seed(SEED)
res <- list(); t0 <- Sys.time(); done <- 0

for (i in ctrl_rows) {
  tr <- draws[[i]]
  sd_i <- diag$seed[i]

  res[[length(res)+1]] <- data.table(
    draw = i, subset_type = "all_36", selector = "-", rep = 0L,
    n_features = d, macro_f1 = fit_eval(tr, seq_len(d), sd_i))

  for (sel in c("L1", "RF")) {
    RM <- if (sel == "L1") sel_obj$rank_l1 else sel_obj$rank_rf
    best  <- RM[i, seq_len(K)]
    worst <- RM[i, (d - K + 1):d]
    res[[length(res)+1]] <- data.table(
      draw = i, subset_type = "selected", selector = sel, rep = 0L,
      n_features = K, macro_f1 = fit_eval(tr, best, sd_i))
    res[[length(res)+1]] <- data.table(
      draw = i, subset_type = "worst", selector = sel, rep = 0L,
      n_features = K, macro_f1 = fit_eval(tr, worst, sd_i))
  }

  for (b in seq_len(N_RANDOM)) {
    rnd <- sample.int(d, K)
    res[[length(res)+1]] <- data.table(
      draw = i, subset_type = "random", selector = "-", rep = b,
      n_features = K, macro_f1 = fit_eval(tr, rnd, sd_i))
  }

  done <- done + 1
  el <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  cat(sprintf("  [%d/%d draws] %.0fs elapsed, ~%.0fs remaining\n",
              done, N_DRAWS, el, el/done * (N_DRAWS - done)))
}

rb <- rbindlist(res)
fwrite(rb, file.path(REP_DIR, "10_random_baseline_raw.csv"))
cat(sprintf("complete in %s\n", format(round(difftime(Sys.time(), t0)))))

cat("\n", strrep("=", 76), "\n", sep = "")
cat("  DOES SELECTION MATTER?  (8-class macro-F1, fixed composition)\n")
cat(strrep("=", 76), "\n", sep = "")

summ <- rb[, .(n = .N, mean = round(mean(macro_f1), 4),
               sd = round(sd(macro_f1), 4),
               min = round(min(macro_f1), 4),
               max = round(max(macro_f1), 4)),
           by = .(subset_type, selector)]
setorder(summ, -mean)
print(summ)
fwrite(summ, file.path(REP_DIR, "10_summary.csv"))

rand_all <- rb[subset_type == "random"]
pct <- rbindlist(lapply(c("L1", "RF"), function(sel) {
  rbindlist(lapply(ctrl_rows, function(i) {
    s <- rb[draw == i & subset_type == "selected" & selector == sel, macro_f1]
    r <- rand_all[draw == i, macro_f1]
    data.table(selector = sel, draw = i, selected = s,
               random_mean = mean(r),
               pct_random_beaten = mean(r < s),
               n_random_better = sum(r >= s))
  }))
}))
cat("\n--- selected subset vs the random distribution, per draw ---\n")
print(pct[, .(selector, selected = round(selected, 4),
              random_mean = round(random_mean, 4),
              gap = round(selected - random_mean, 4),
              pct_random_beaten = round(100 * pct_random_beaten, 1))], nrows = 30)

cat("\noverall \n")
ov <- pct[, .(draws = .N,
              mean_selected = round(mean(selected), 4),
              mean_random   = round(mean(random_mean), 4),
              mean_gap      = round(mean(selected - random_mean), 4),
              mean_pct_random_beaten = round(100 * mean(pct_random_beaten), 1)),
          by = selector]
print(ov)
fwrite(pct, file.path(REP_DIR, "10_selected_vs_random.csv"))

for (sel in c("L1", "RF")) {
  p <- pct[selector == sel]
  w <- suppressWarnings(wilcox.test(p$selected, p$random_mean,
                                    paired = TRUE, alternative = "greater"))
  cat(sprintf("  %s: median gap %+.4f, V = %s, p = %.4f\n",
              sel, median(p$selected - p$random_mean),
              w$statistic, w$p.value))
}

gap_L1 <- ov[selector == "L1", mean_gap]
gap_RF <- ov[selector == "RF", mean_gap]
rng <- summ[subset_type == "all_36", mean] - summ[subset_type == "worst", mean][1]

cat("\n", strrep("-", 76), "\n", sep = "")
cat(sprintf("  all-36 mean            : %.4f\n", summ[subset_type == "all_36", mean]))
cat(sprintf("  selected mean (L1/RF)  : %.4f / %.4f\n",
            summ[subset_type == "selected" & selector == "L1", mean],
            summ[subset_type == "selected" & selector == "RF", mean]))
cat(sprintf("  random mean            : %.4f\n", summ[subset_type == "random", mean]))
cat(sprintf("  worst mean (L1/RF)     : %.4f / %.4f\n",
            summ[subset_type == "worst" & selector == "L1", mean],
            summ[subset_type == "worst" & selector == "RF", mean]))
cat(strrep("-", 76), "\n", sep = "")

cat("\nINTERPRETATION\n")
if (max(gap_L1, gap_RF) < 0.01) {
  cat("  selected ~= random. The feature space is REDUNDANT: many different\n")
  cat("  10-subsets carry the same information. Report this plainly. The\n")
  cat("  Stage 8 stability finding still stands - it concerns WHICH features\n")
  cat("  are named as informative (interpretation, and what a deployed sensor\n")
  cat("  must collect), not predictive performance.\n")
} else if (max(gap_L1, gap_RF) < 0.05) {
  cat("  selected > random, but modestly. Selection helps; which good subset\n")
  cat("  you land on matters much less. This explains the near-zero Stage 9\n")
  cat("  cost without appealing to redundancy.\n")
} else {
  cat("  selected >> random. Selection matters substantially, so the near-zero\n")
  cat("  Stage 9 cost needs another explanation - most likely that the\n")
  cat("  composition perturbations never pushed the subset out of the good\n")
  cat("  region of the feature space.\n")
}


## FIGURES 33-34

lab <- c(all_36 = "All 36 features", selected = "Selected top-10",
         random = "Random 10", worst = "Worst 10 (bottom-ranked)")
rb[, type_lab := factor(lab[subset_type], levels = rev(lab))]

p33 <- ggplot(rb, aes(type_lab, macro_f1, fill = subset_type)) +
  geom_boxplot(alpha = 0.85, outlier.size = 0.7, show.legend = FALSE) +
  geom_jitter(width = 0.12, alpha = 0.3, size = 0.9, show.legend = FALSE) +
  coord_flip() +
  scale_fill_manual(values = c(all_36 = "#2E7D32", selected = "#1565C0",
                               random = "grey60", worst = "#C62828")) +
  labs(title = "Fig 33. Does feature selection matter on this dataset?",
       subtitle = sprintf("8-class macro-F1, %d control draws at fixed composition, same model throughout.\nIf 'Selected' sits on top of 'Random', many different subsets carry the same information.",
                          N_DRAWS),
       x = NULL, y = "macro-F1 on the held-out test partition") +
  theme_ids()
save_fig(p33, "33_selection_vs_random", w = 10, h = 5)

p34 <- ggplot(rand_all, aes(macro_f1)) +
  geom_histogram(bins = 25, fill = "grey65", colour = "white") +
  geom_vline(data = rb[subset_type == "selected"],
             aes(xintercept = macro_f1, colour = selector), linewidth = 0.6) +
  geom_vline(xintercept = mean(rb[subset_type == "all_36", macro_f1]),
             colour = "#2E7D32", linewidth = 0.9, linetype = "dashed") +
  scale_colour_manual(values = c(L1 = "#1565C0", RF = "#EF6C00")) +
  labs(title = sprintf("Fig 34. The selected subsets against %d random 10-feature subsets",
                       nrow(rand_all)),
       subtitle = "Grey: random subsets. Coloured lines: the selectors' choices. Green dashed: all 36 features.\nHow far right the coloured lines sit is how much selection actually buys.",
       x = "macro-F1", y = "random subsets", colour = "Selector") +
  theme_ids()
save_fig(p34, "34_random_distribution", w = 10, h = 5)

saveRDS(list(raw = rb, summary = summ, per_draw = pct, overall = ov,
             k = K, n_draws = N_DRAWS, n_random = N_RANDOM,
             num_trees = NUM_TREES, test_n = length(te_idx), seed = SEED),
        file.path(OUT_DIR, "processed", "random_baseline.rds"))

si <- c(
  sprintf("timestamp     : %s", format(Sys.time())),
  sprintf("R version     : %s", R.version.string),
  sprintf("seed          : %d", SEED),
  sprintf("k             : %d", K),
  sprintf("control draws : %d", N_DRAWS),
  sprintf("random subsets: %d per draw (%d total)", N_RANDOM, nrow(rand_all)),
  sprintf("test rows     : %d", length(te_idx)),
  sprintf("all_36 mean   : %.4f", summ[subset_type == "all_36", mean]),
  sprintf("selected L1   : %.4f", summ[subset_type == "selected" & selector == "L1", mean]),
  sprintf("selected RF   : %.4f", summ[subset_type == "selected" & selector == "RF", mean]),
  sprintf("random mean   : %.4f", summ[subset_type == "random", mean]),
  sprintf("worst L1      : %.4f", summ[subset_type == "worst" & selector == "L1", mean]),
  sprintf("worst RF      : %.4f", summ[subset_type == "worst" & selector == "RF", mean])
)
writeLines(si, file.path(REP_DIR, "10_run_record.txt"))
cat("\n", paste(si, collapse = "\n"), "\n")
cat("\nStage 10 complete. -> data/processed/random_baseline.rds\n")
