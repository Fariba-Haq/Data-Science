

SEED        <- 13
NUM_TREES   <- 100        # ranger trees per evaluation model
TEST_SAMPLE <- 60000      # stratified test rows used for evaluation
ARMS_TO_RUN <- c("control", "drop_one", "dirichlet", "reference")

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
arms    <- readRDS(file.path(OUT_DIR, "processed", "arms.rds"))
sel_obj <- readRDS(file.path(OUT_DIR, "processed", "selection.rds"))

Xtr <- sp$X_train; y8tr <- sp$y8_train
Xte <- sp$X_test;  y8te <- sp$y8_test
FEATS <- sp$features; d <- length(FEATS)
K_MAIN <- sel_obj$k_main
safe_names <- make.names(FEATS, unique = TRUE)

diag  <- sel_obj$diag
keep  <- diag$arm %in% ARMS_TO_RUN
rows_use <- which(keep)
draws <- arms$draws[diag$draw_key]

set.seed(SEED)
te_idx <- unlist(lapply(unique(y8te), function(cl) {
  w <- which(y8te == cl)
  sample(w, max(1, round(length(w) * TEST_SAMPLE / length(y8te))))
}))
Xte_s <- Xte[te_idx, , drop = FALSE]; yte_s <- factor(y8te[te_idx])
cat(sprintf("Stage 9: %d configurations | test subsample %s rows\n",
            length(rows_use), format(length(te_idx), big.mark = ",")))
print(table(yte_s))

macro_f1 <- function(truth, pred, levs) {
  truth <- factor(truth, levels = levs); pred <- factor(pred, levels = levs)
  f1 <- sapply(levs, function(cl) {
    tp <- sum(pred == cl & truth == cl)
    fp <- sum(pred == cl & truth != cl)
    fn <- sum(pred != cl & truth == cl)
    if (tp == 0) return(0)
    p <- tp / (tp + fp); r <- tp / (tp + fn)
    if (p + r == 0) 0 else 2 * p * r / (p + r)
  })
  list(macro = mean(f1), per_class = f1)
}

LEVS <- levels(yte_s)

fit_eval <- function(train_rows, cols, seed) {
  x <- Xtr[train_rows, cols, drop = FALSE]
  colnames(x) <- safe_names[cols]
  rf <- ranger(x = x, y = factor(y8tr[train_rows], levels = LEVS),
               num.trees = NUM_TREES, num.threads = 0,
               seed = seed, verbose = FALSE)
  xt <- Xte_s[, cols, drop = FALSE]; colnames(xt) <- safe_names[cols]
  pr <- predict(rf, data = xt, num.threads = 0)$predictions
  macro_f1(yte_s, pr, LEVS)
}

consensus <- function(rank_mat, rws, k) {
  freq <- integer(d); rsum <- numeric(d)
  for (r in rws) {
    s <- rank_mat[r, seq_len(k)]
    freq[s] <- freq[s] + 1L
    rsum <- rsum + order(rank_mat[r, ])
  }
  order(-freq, rsum)[seq_len(k)]
}
ctrl_rows <- which(diag$config == "control_base")
ctrl_set <- list(
  L1 = consensus(sel_obj$rank_l1, ctrl_rows, K_MAIN),
  RF = consensus(sel_obj$rank_rf, ctrl_rows, K_MAIN))
cat("\ncontrol consensus (committed) subsets:\n")
for (s in names(ctrl_set))
  cat(sprintf("  %s: %s\n", s, paste(FEATS[sort(ctrl_set[[s]])], collapse = ", ")))

res <- list(); t0 <- Sys.time(); done <- 0
for (i in rows_use) {
  tr <- draws[[i]]
  a_all <- fit_eval(tr, seq_len(d), diag$seed[i])          # 36-feature reference

  for (sel in c("L1", "RF")) {
    RM  <- if (sel == "L1") sel_obj$rank_l1 else sel_obj$rank_rf
    own <- RM[i, seq_len(K_MAIN)]
    ctl <- ctrl_set[[sel]]

    a_own <- fit_eval(tr, own, diag$seed[i])
    a_ctl <- if (setequal(own, ctl)) a_own else fit_eval(tr, ctl, diag$seed[i])

    res[[length(res) + 1]] <- data.table(
      selector = sel, draw_key = diag$draw_key[i], arm = diag$arm[i],
      config = diag$config[i], js = diag$js_realised_vs_base[i],
      f1_own = a_own$macro, f1_ctrl = a_ctl$macro, f1_all = a_all$macro,
      selection_cost = a_own$macro - a_ctl$macro,
      reduction_cost = a_all$macro - a_own$macro,
      overlap_own_ctrl = length(intersect(own, ctl)) /
                         length(union(own, ctl)))
  }
  done <- done + 1
  if (done %% 10 == 0 || done == length(rows_use)) {
    el <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    cat(sprintf("  [%d/%d] %.0fs elapsed, ~%.0fs remaining\n",
                done, length(rows_use), el, el / done * (length(rows_use) - done)))
  }
}
f1 <- rbindlist(res)
fwrite(f1, file.path(REP_DIR, "09_macro_f1_per_draw.csv"))
cat(sprintf("complete in %s\n", format(round(difftime(Sys.time(), t0)))))


cat("\n", strrep("=", 78), "\n", sep = "")
cat("  MACRO-F1 BY ARM  (8-class, held-out test partition)\n")
cat(strrep("=", 78), "\n", sep = "")
by_arm <- f1[, .(draws = .N,
                 f1_all  = round(mean(f1_all), 4),
                 f1_own  = round(mean(f1_own), 4),
                 f1_ctrl = round(mean(f1_ctrl), 4),
                 selection_cost = round(mean(selection_cost), 4),
                 cost_sd = round(sd(selection_cost), 4)),
             by = .(selector, arm)]
setorder(by_arm, selector, -f1_ctrl)
print(by_arm)
fwrite(by_arm, file.path(REP_DIR, "09_macro_f1_by_arm.csv"))

cat("\n--- reduction cost: 36 features -> top", K_MAIN, "---\n")
print(f1[, .(f1_all = round(mean(f1_all), 4),
             f1_own = round(mean(f1_own), 4),
             reduction_cost = round(mean(reduction_cost), 4)), by = selector])

cat("\n--- selection cost vs composition shift ---\n")
cc <- f1[js > 0, .(n = .N,
                   r_cost_js = round(cor(js, selection_cost), 3),
                   r_cost_overlap = round(cor(overlap_own_ctrl, selection_cost), 3)),
         by = selector]
print(cc)
cat("\nr_cost_overlap should be NEGATIVE: the less the committed subset\n")
cat("overlaps the one the mix would have chosen, the more macro-F1 is lost.\n")

cat("\n--- drop-one, ranked by what committing to the control subset costs ---\n")
print(f1[arm %in% c("drop_one", "reference"),
         .(f1_own = round(mean(f1_own), 4),
           f1_ctrl = round(mean(f1_ctrl), 4),
           cost = round(mean(selection_cost), 4)),
         by = .(selector, config)][order(selector, -cost)], nrows = 20)

p30 <- ggplot(f1, aes(arm, f1_ctrl, fill = arm)) +
  geom_boxplot(alpha = 0.85, outlier.size = 0.6, show.legend = FALSE) +
  facet_wrap(~selector) + coord_flip() +
  labs(title = sprintf("Fig 30. Macro-F1 of the committed subset as the threat mix moves (k = %d)", K_MAIN),
       subtitle = "The detector commits to the control-consensus subset, then the attack mix changes.\n8-class macro-F1 on the held-out test partition.",
       x = NULL, y = "macro-F1 using the committed subset") +
  theme_ids()
save_fig(p30, "30_macro_f1_by_arm", w = 10, h = 5)

p31 <- ggplot(f1[js > 0], aes(js, selection_cost, colour = arm)) +
  geom_hline(yintercept = 0, colour = "grey50") +
  geom_point(alpha = 0.5, size = 1.4) +
  geom_smooth(method = "loess", formula = y ~ x, se = TRUE,
              colour = "#C62828", linewidth = 0.8) +
  facet_wrap(~selector) +
  labs(title = "Fig 31. The macro-F1 cost of a committed feature subset",
       subtitle = "Cost = macro-F1 with the mix's own subset minus macro-F1 with the committed subset.\nSame rows, same model, same k - only the chosen features differ, so the gap is attributable to selection.",
       x = "Jensen-Shannon divergence from the base mix",
       y = "selection cost (macro-F1)", colour = NULL) +
  theme_ids()
save_fig(p31, "31_selection_cost", w = 10, h = 5)

p32 <- ggplot(f1[js > 0], aes(overlap_own_ctrl, selection_cost, colour = arm)) +
  geom_hline(yintercept = 0, colour = "grey50") +
  geom_point(alpha = 0.5, size = 1.4) +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE,
              colour = "#1565C0", linewidth = 0.8) +
  facet_wrap(~selector) +
  labs(title = "Fig 32. Subset overlap predicts the macro-F1 cost",
       subtitle = "The link between the stability result and the performance result: subsets that diverge\nfrom the committed one cost more macro-F1. This is why subset overlap is worth measuring.",
       x = "Jaccard between the mix's own subset and the committed subset",
       y = "selection cost (macro-F1)", colour = NULL) +
  theme_ids()
save_fig(p32, "32_overlap_vs_cost", w = 10, h = 5)


saveRDS(list(f1 = f1, by_arm = by_arm, ctrl_set = ctrl_set,
             k = K_MAIN, num_trees = NUM_TREES, test_n = length(te_idx),
             seed = SEED),
        file.path(OUT_DIR, "processed", "macro_f1.rds"))

si <- c(
  sprintf("timestamp     : %s", format(Sys.time())),
  sprintf("R version     : %s", R.version.string),
  sprintf("ranger        : %s", as.character(packageVersion("ranger"))),
  sprintf("seed          : %d", SEED),
  sprintf("k             : %d", K_MAIN),
  sprintf("trees         : %d", NUM_TREES),
  sprintf("test rows     : %d", length(te_idx)),
  sprintf("configurations: %d", length(rows_use)),
  sprintf("mean f1_all   : %.4f", mean(f1$f1_all)),
  sprintf("mean f1_own   : %.4f", mean(f1$f1_own)),
  sprintf("mean f1_ctrl  : %.4f", mean(f1$f1_ctrl)),
  sprintf("mean selection cost: %.4f", mean(f1$selection_cost))
)
writeLines(si, file.path(REP_DIR, "09_run_record.txt"))
cat("\n", paste(si, collapse = "\n"), "\n")
cat("\nStage 9 complete. -> data/processed/macro_f1.rds\n")
