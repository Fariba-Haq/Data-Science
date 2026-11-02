
SEED       <- 13
N_DRAWS    <- 4          # control draws (composition held fixed)
N_REPS     <- 5          # random subsets per overlap level per draw
SELECTORS  <- c("L1")    # add "RF" to run both (doubles the runtime)
NUM_TREES  <- 100
TEST_SAMPLE<- 40000

OUT_DIR <- "data"; FIG_DIR <- "figures"; REP_DIR <- "reports"

need <- c("data.table", "ggplot2", "scales", "ranger")
miss <- need[!sapply(need, requireNamespace, quietly = TRUE)]
if (length(miss))
  stop("Install first: install.packages(c(",
       paste(sprintf("'%s'", miss), collapse = ", "), "))")
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(scales); library(ranger)
})
for (d_ in c(FIG_DIR, REP_DIR)) dir.create(d_, showWarnings = FALSE, recursive = TRUE)

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
rbase   <- tryCatch(readRDS(file.path(OUT_DIR, "processed", "random_baseline.rds")),
                    error = function(e) NULL)

Xtr <- sp$X_train; y8tr <- sp$y8_train
Xte <- sp$X_test;  y8te <- sp$y8_test
FEATS <- sp$features; d <- length(FEATS)
K <- sel_obj$k_main
safe_names <- make.names(FEATS, unique = TRUE)

diag  <- sel_obj$diag
draws <- arms$draws[diag$draw_key]
ctrl_rows_all <- which(diag$config == "control_base")
use_draws <- ctrl_rows_all[seq_len(N_DRAWS)]

set.seed(SEED)
te_idx <- unlist(lapply(unique(y8te), function(cl) {
  w <- which(y8te == cl)
  sample(w, max(1, round(length(w) * TEST_SAMPLE / length(y8te))))
}))
Xte_s <- Xte[te_idx, , drop = FALSE]; yte_s <- factor(y8te[te_idx])
LEVS <- levels(yte_s)

n_fits <- N_DRAWS * (1 + length(SELECTORS) * (K + 1) * N_REPS)
cat(sprintf("Stage 11: %d draws x %d selector(s) x %d overlap levels x %d reps\n",
            N_DRAWS, length(SELECTORS), K + 1, N_REPS))
cat(sprintf("          ~%d model fits | test %s rows | k = %d\n",
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

consensus <- function(rank_mat, rws, k) {
  freq <- integer(d); rsum <- numeric(d)
  for (r in rws) {
    s <- rank_mat[r, seq_len(k)]
    freq[s] <- freq[s] + 1L
    rsum <- rsum + order(rank_mat[r, ])
  }
  order(-freq, rsum)[seq_len(k)]
}
ctrl_set <- list(
  L1 = consensus(sel_obj$rank_l1, ctrl_rows_all, K),
  RF = consensus(sel_obj$rank_rf, ctrl_rows_all, K))

cat("\ncommitted subsets:\n")
for (s in SELECTORS)
  cat(sprintf("  %s: %s\n", s, paste(FEATS[sort(ctrl_set[[s]])], collapse = ", ")))

set.seed(SEED)
res <- list(); t0 <- Sys.time(); done <- 0

for (i in use_draws) {
  tr <- draws[[i]]; sd_i <- diag$seed[i]

  res[[length(res)+1]] <- data.table(
    draw = i, selector = "-", m = NA_integer_, rep = 0L,
    jaccard = NA_real_, subset_type = "all_36",
    features = paste(sort(FEATS), collapse = "|"),
    macro_f1 = fit_eval(tr, seq_len(d), sd_i))

  for (sel in SELECTORS) {
    C <- ctrl_set[[sel]]
    outside <- setdiff(seq_len(d), C)

    for (m in 0:K) {
      for (r in seq_len(N_REPS)) {
       
        keep_in  <- if (m > 0) sample(C, m) else integer(0)
        keep_out <- if (K - m > 0) sample(outside, K - m) else integer(0)
        cols <- c(keep_in, keep_out)

        res[[length(res)+1]] <- data.table(
          draw = i, selector = sel, m = m, rep = r,
          jaccard = m / (2*K - m),
          subset_type = if (m == K) "committed" else "constructed",
          features = paste(sort(FEATS[cols]), collapse = "|"),
          macro_f1 = fit_eval(tr, cols, sd_i))

       
        if (m == K) break
      }
    }
  }

  done <- done + 1
  el <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  cat(sprintf("  [%d/%d draws] %.0fs elapsed, ~%.0fs remaining\n",
              done, N_DRAWS, el, el/done * (N_DRAWS - done)))
}

oc <- rbindlist(res)
fwrite(oc, file.path(REP_DIR, "11_overlap_curve_raw.csv"))
cat(sprintf("complete in %s\n", format(round(difftime(Sys.time(), t0)))))

cur <- oc[!is.na(m), .(n = .N,
                       mean = round(mean(macro_f1), 4),
                       sd   = round(sd(macro_f1), 4),
                       lo   = round(quantile(macro_f1, 0.025, names = FALSE), 4),
                       hi   = round(quantile(macro_f1, 0.975, names = FALSE), 4)),
          by = .(selector, m, jaccard)]
setorder(cur, selector, -m)

cat("\n", strrep("=", 74), "\n", sep = "")
cat("  MACRO-F1 AGAINST OVERLAP WITH THE COMMITTED SUBSET\n")
cat(strrep("=", 74), "\n", sep = "")
print(cur[, .(selector, features_kept = m, jaccard = round(jaccard, 3),
              n, mean, sd, ci = sprintf("[%.3f, %.3f]", lo, hi))], nrows = 30)
fwrite(cur, file.path(REP_DIR, "11_overlap_curve.csv"))

TOL <- 0.005
cat("\n--- where performance starts to fall (tolerance", TOL, ") ---\n")
for (sel in SELECTORS) {
  cs <- cur[selector == sel][order(-m)]
  ref <- cs[m == K, mean]
  drop <- cs[mean < ref - TOL]
  first <- if (nrow(drop)) drop[1] else NULL
  if (is.null(first)) {
    cat(sprintf("  %s: performance never falls more than %.3f - flat throughout\n",
                sel, TOL))
  } else {
    cat(sprintf("  %s: committed = %.4f | first drop at m = %d (Jaccard %.3f), mean %.4f\n",
                sel, ref, first$m, first$jaccard, first$mean))
    cat(sprintf("      -> the plateau runs from Jaccard 1.000 down to %.3f\n",
                cs[m == first$m + 1, jaccard]))
  }
}

cat("\n--- Stage 8 observed this range of overlaps ---\n")
cat("  control vs control      : 0.86\n")
cat("  composition-perturbed   : 0.59 - 0.86\n")
cat("  random (Stage 10, expected): ~0.16\n")
cat("Compare those against the curve above: if 0.59-0.86 sits on the plateau\n")
cat("and ~0.16 sits past the cliff, Stages 8, 9 and 10 are one coherent result.\n")


ref_all  <- mean(oc[subset_type == "all_36", macro_f1])
ref_rand <- if (!is.null(rbase)) rbase$summary[subset_type == "random", mean] else NA_real_

p35 <- ggplot(oc[!is.na(m)], aes(jaccard, macro_f1)) +
  geom_hline(yintercept = ref_all, colour = "#2E7D32",
             linetype = "dashed", linewidth = 0.7) +
  { if (!is.na(ref_rand))
      geom_hline(yintercept = ref_rand, colour = "grey45", linetype = "dotted",
                 linewidth = 0.7) } +
  annotate("rect", xmin = 0.59, xmax = 0.86, ymin = -Inf, ymax = Inf,
           fill = "#1565C0", alpha = 0.10) +
  geom_point(alpha = 0.45, size = 1.5, colour = "#1565C0") +
  geom_smooth(method = "loess", formula = y ~ x, se = TRUE,
              colour = "#C62828", linewidth = 0.9) +
  facet_wrap(~selector) +
  scale_x_continuous(breaks = seq(0, 1, 0.2)) +
  labs(title = "Fig 35. The plateau and the cliff",
       subtitle = paste0("Each point is a subset built to have EXACTLY m of the committed subset's 10 features.\n",
                         "Shaded band: the overlap range the composition arms actually produced (Stage 8).\n",
                         "Green dashed: all 36 features. Grey dotted: random 10-subsets (Stage 10)."),
       x = "Jaccard overlap with the committed subset",
       y = "8-class macro-F1") +
  theme_ids()
save_fig(p35, "35_plateau_and_cliff", w = 10, h = 5.5)

p36 <- ggplot(cur, aes(factor(m), mean)) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.25, colour = "grey40") +
  geom_point(size = 2.6, colour = "#1565C0") +
  geom_hline(yintercept = ref_all, colour = "#2E7D32", linetype = "dashed") +
  facet_wrap(~selector) +
  labs(title = "Fig 36. Macro-F1 by number of committed features retained",
       subtitle = paste0("m = how many of the committed subset's 10 features the model was given;\n",
                         "the remaining 10 - m were drawn from the other 26. Green dashed: all 36 features."),
       x = "features retained from the committed subset (m)",
       y = "macro-F1 (mean, 95% interval)") +
  theme_ids()
save_fig(p36, "36_overlap_levels", w = 10, h = 5)

saveRDS(list(raw = oc, curve = cur, ctrl_set = ctrl_set, k = K,
             n_draws = N_DRAWS, n_reps = N_REPS, selectors = SELECTORS,
             num_trees = NUM_TREES, test_n = length(te_idx), seed = SEED),
        file.path(OUT_DIR, "processed", "overlap_curve.rds"))

si <- c(
  sprintf("timestamp    : %s", format(Sys.time())),
  sprintf("R version    : %s", R.version.string),
  sprintf("seed         : %d", SEED),
  sprintf("k            : %d", K),
  sprintf("draws        : %d", N_DRAWS),
  sprintf("reps / level : %d", N_REPS),
  sprintf("selectors    : %s", paste(SELECTORS, collapse = ", ")),
  sprintf("test rows    : %d", length(te_idx)),
  sprintf("all_36 mean  : %.4f", ref_all),
  sprintf("committed    : %s",
          paste(sprintf("%s=%.4f", SELECTORS,
                        sapply(SELECTORS, function(s) cur[selector == s & m == K, mean])),
                collapse = " ")),
  sprintf("m=0 (no overlap): %s",
          paste(sprintf("%s=%.4f", SELECTORS,
                        sapply(SELECTORS, function(s) cur[selector == s & m == 0, mean])),
                collapse = " "))
)
writeLines(si, file.path(REP_DIR, "11_run_record.txt"))
cat("\n", paste(si, collapse = "\n"), "\n")
cat("\nStage 11 complete. -> data/processed/overlap_curve.rds\n")
