

SEED        <- 13
K_VALUES    <- c(5, 10, 15)     
K_MAIN      <- 10               
NUM_TREES   <- 200              
ARMS_TO_RUN <- c("control", "drop_one", "dirichlet", "reference")
                                

OUT_DIR <- "data"; FIG_DIR <- "figures"; REP_DIR <- "reports"

need <- c("data.table", "ggplot2", "scales", "glmnet", "ranger")
miss <- need[!sapply(need, requireNamespace, quietly = TRUE)]
if (length(miss))
  stop("Install first:  install.packages(c(",
       paste(sprintf("'%s'", miss), collapse = ", "), "))")

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(scales)
  library(glmnet); library(ranger)
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


sp   <- readRDS(file.path(OUT_DIR, "processed", "split_safe.rds"))
arms <- readRDS(file.path(OUT_DIR, "processed", "arms.rds"))

X    <- sp$X_train
ybin <- sp$ybin_train
FEATS <- sp$features
d     <- length(FEATS)

safe_names <- make.names(FEATS, unique = TRUE)
name_map   <- setNames(FEATS, safe_names)

diag <- arms$diag
keep <- diag$arm %in% ARMS_TO_RUN
diag <- diag[keep]
draws <- arms$draws[diag$draw_key]

cat(sprintf("Stage 7: %d configurations x 2 selectors, %d features, k = %s\n",
            nrow(diag), d, paste(K_VALUES, collapse = ", ")))
print(diag[, .(configs = .N), by = arm])

select_l1 <- function(x, y) {
  fit <- glmnet(x, y, family = "binomial", alpha = 1,
                standardize = FALSE, nlambda = 60)
  B <- as.matrix(fit$beta) != 0
  entry <- apply(B, 1, function(r) { w <- which(r)[1]; if (length(w)) w else Inf })
  final <- abs(as.numeric(fit$beta[, ncol(fit$beta)]))

  order(entry, -final)
}

select_rf <- function(x, y, seed) {
  colnames(x) <- safe_names
  rf <- ranger(x = x, y = factor(y), num.trees = NUM_TREES,
               importance = "impurity", num.threads = 0,
               seed = seed, verbose = FALSE)
  order(rf$variable.importance, decreasing = TRUE)
}


rank_l1 <- matrix(NA_integer_, nrow = nrow(diag), ncol = d)
rank_rf <- matrix(NA_integer_, nrow = nrow(diag), ncol = d)
t0 <- Sys.time()

for (i in seq_len(nrow(diag))) {
  idx <- draws[[i]]
  xi  <- X[idx, , drop = FALSE]
  yi  <- ybin[idx]

  rank_l1[i, ] <- select_l1(xi, yi)
  rank_rf[i, ] <- select_rf(xi, yi, seed = diag$seed[i])

  rm(xi, yi)
  if (i %% 10 == 0 || i == nrow(diag)) {
    el  <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    eta <- el / i * (nrow(diag) - i)
    cat(sprintf("  [%d/%d] %.0fs elapsed, ~%.0fs remaining\n",
                i, nrow(diag), el, eta))
  }
}
cat(sprintf("selection complete in %s\n",
            format(round(difftime(Sys.time(), t0)))))

topk <- function(rank_mat, rows, k) {
  Z <- matrix(0L, nrow = length(rows), ncol = d)
  for (j in seq_along(rows)) Z[j, rank_mat[rows[j], seq_len(k)]] <- 1L
  Z
}


nogueira <- function(Z) {
  M <- nrow(Z); dd <- ncol(Z)
  if (M < 2) return(NA_real_)
  pf   <- colMeans(Z)
  kbar <- mean(rowSums(Z))
  num  <- mean((M / (M - 1)) * pf * (1 - pf))
  den  <- (kbar / dd) * (1 - kbar / dd)
  if (den == 0) return(NA_real_)
  1 - num / den
}

jaccard <- function(a, b) {
  u <- length(union(a, b)); if (u == 0) return(NA_real_)
  length(intersect(a, b)) / u
}

consensus <- function(rank_mat, rows, k) {
  freq <- integer(d); rsum <- numeric(d)
  for (r in rows) {
    sel <- rank_mat[r, seq_len(k)]
    freq[sel] <- freq[sel] + 1L
    rsum <- rsum + order(rank_mat[r, ])
  }
  order(-freq, rsum)[seq_len(k)]
}

res_stab <- list(); res_jac <- list()
configs  <- unique(diag$config)

for (sel in c("L1", "RF")) {
  RM <- if (sel == "L1") rank_l1 else rank_rf
  for (k in K_VALUES) {
    ctrl_rows <- which(diag$config == "control_base")
    ctrl_cons <- if (length(ctrl_rows)) consensus(RM, ctrl_rows, k) else NULL

    for (cf in configs) {
      rows <- which(diag$config == cf)
      Z    <- topk(RM, rows, k)
      cons <- consensus(RM, rows, k)

      res_stab[[length(res_stab) + 1]] <- data.table(
        selector = sel, k = k, config = cf,
        arm = diag$arm[rows[1]], n_draws = length(rows),
        js  = diag$js_realised_vs_base[rows[1]],
        nogueira = nogueira(Z))

      if (!is.null(ctrl_cons))
        res_jac[[length(res_jac) + 1]] <- data.table(
          selector = sel, k = k, config = cf,
          arm = diag$arm[rows[1]],
          js  = diag$js_realised_vs_base[rows[1]],
          jaccard_vs_control = jaccard(cons, ctrl_cons),
          consensus = paste(FEATS[sort(cons)], collapse = "|"))
    }
  }
}
stab <- rbindlist(res_stab)
jac  <- rbindlist(res_jac)

fwrite(stab, file.path(REP_DIR, "07_stability.csv"))
fwrite(jac,  file.path(REP_DIR, "07_jaccard_vs_control.csv"))


cat("\n", strrep("=", 74), "\n", sep = "")
cat("  CONTROL ARM — THE FALSIFIER\n")
cat(strrep("=", 74), "\n", sep = "")
ctrl <- stab[config == "control_base"]
print(ctrl[, .(selector, k, n_draws, nogueira = round(nogueira, 4))])

for (sel in unique(ctrl$selector)) {
  v <- ctrl[selector == sel & k == K_MAIN, nogueira]
  if (!length(v) || is.na(v)) next
  verdict <- if (v >= 0.8) "STABLE - design holds, movement elsewhere is attributable to composition"
        else if (v >= 0.6) "MARGINAL - report the noise floor prominently and interpret cautiously"
        else               "UNSTABLE - selection moves at FIXED composition. STOP and reconsider."
  cat(sprintf("\n  %s at k=%d : %.4f  -> %s\n", sel, K_MAIN, v, verdict))
}

if (nrow(jac)) {
  cat("\n", strrep("-", 74), "\n", sep = "")
  cat("  SUBSET OVERLAP vs CONTROL CONSENSUS  (k =", K_MAIN, ")\n")
  cat(strrep("-", 74), "\n", sep = "")
  print(jac[k == K_MAIN & arm %in% c("control", "drop_one", "reference")][
        order(selector, -jaccard_vs_control),
        .(selector, config, js = round(js, 4),
          jaccard = round(jaccard_vs_control, 4))], nrows = 40)

  cat("\n--- per-arm means ---\n")
  print(jac[k == K_MAIN, .(configs = .N,
                           jaccard_mean = round(mean(jaccard_vs_control), 4),
                           jaccard_min  = round(min(jaccard_vs_control), 4)),
            by = .(selector, arm)])
}

cat("\n--- control consensus subsets (k =", K_MAIN, ") ---\n")
print(jac[config == "control_base" & k == K_MAIN, .(selector, consensus)])

cs <- jac[config == "control_base"]
if (uniqueN(cs$selector) == 2) {
  cat("\n--- L1 vs RF agreement at fixed composition ---\n")
  agree <- rbindlist(lapply(K_VALUES, function(kk) {
    a <- strsplit(cs[selector == "L1" & k == kk, consensus], "|", fixed = TRUE)[[1]]
    b <- strsplit(cs[selector == "RF" & k == kk, consensus], "|", fixed = TRUE)[[1]]
    data.table(k = kk, shared = length(intersect(a, b)),
               jaccard_L1_vs_RF = round(jaccard(a, b), 4),
               only_L1 = paste(setdiff(a, b), collapse = ", "),
               only_RF = paste(setdiff(b, a), collapse = ", "))
  }))
  print(agree)
  fwrite(agree, file.path(REP_DIR, "07_selector_agreement.csv"))
}


## Fig 22: the falsifier, in one picture
p22 <- ggplot(stab[!is.na(nogueira)],
              aes(factor(k), nogueira, fill = arm)) +
  geom_boxplot(outlier.size = 0.6, alpha = 0.85) +
  facet_wrap(~selector) +
  geom_hline(yintercept = 0.8, linetype = "dashed", colour = "grey40") +
  coord_cartesian(ylim = c(0, 1)) +
  labs(title = "Fig 22. Within-configuration selection stability (Nogueira index)",
       subtitle = "The control arm is the noise floor: same composition, different draw. If the other arms sit\nbelow it, composition is moving the subset. Dashed line marks 0.8.",
       x = "top-k budget", y = "Nogueira stability", fill = NULL) +
  theme_ids()
save_fig(p22, "22_stability_by_arm", w = 10, h = 5.5)

## Fig 23: the headline arm
if (nrow(jac[k == K_MAIN & arm %in% c("drop_one", "control", "reference")])) {
  dj <- jac[k == K_MAIN & arm %in% c("drop_one", "control", "reference")]
  p23 <- ggplot(dj, aes(reorder(config, jaccard_vs_control),
                        jaccard_vs_control, fill = arm)) +
    geom_col(position = "dodge") + coord_flip() +
    facet_wrap(~selector) +
    scale_y_continuous(limits = c(0, 1)) +
    labs(title = sprintf("Fig 23. Subset overlap with the control consensus (k = %d)", K_MAIN),
         subtitle = "Jaccard of each configuration's consensus subset against the control consensus.\nThe control bar is 1.0 by definition; every drop-one perturbation is of identical magnitude (JS = 0.075).",
         x = NULL, y = "Jaccard vs control", fill = NULL) +
    theme_ids()
  save_fig(p23, "23_jaccard_drop_one", w = 11, h = 6)

  ##  Fig 24: dose-response.
  dd <- jac[k == K_MAIN & arm == "dirichlet"]
  if (nrow(dd)) {
    p24 <- ggplot(dd, aes(js, jaccard_vs_control)) +
      geom_point(alpha = 0.6, size = 1.6, colour = "#6A1B9A") +
      geom_smooth(method = "loess", formula = y ~ x, se = TRUE,
                  colour = "#C62828", linewidth = 0.8) +
      facet_wrap(~selector) +
      scale_y_continuous(limits = c(0, 1)) +
      labs(title = sprintf("Fig 24. Dose-response: subset overlap against composition shift (k = %d)", K_MAIN),
           subtitle = "Pre-registered prediction: the relationship is hyperbolic rather than linear\n(Lukaszuk & Krawczuk 2024; Pes 2020).",
           x = "Jensen-Shannon divergence from the base mix",
           y = "Jaccard vs control consensus") +
      theme_ids()
    save_fig(p24, "24_dose_response", w = 10, h = 5)
  } else {
    cat("  (Fig 24 skipped: no dirichlet configurations in this run)\n")
  }
}

## Fig 25: which features are actually stable
freq_dt <- rbindlist(lapply(c("L1", "RF"), function(sel) {
  RM <- if (sel == "L1") rank_l1 else rank_rf
  cr <- which(diag$config == "control_base")
  ot <- which(diag$arm %in% c("drop_one", "dirichlet"))
  f  <- function(rows) if (!length(rows)) rep(NA_real_, d) else
          colMeans(topk(RM, rows, K_MAIN))
  data.table(selector = sel, feature = FEATS,
             control = f(cr), perturbed = f(ot))
}))
fl <- melt(freq_dt, id.vars = c("selector", "feature"),
           variable.name = "arm", value.name = "freq")
fl <- fl[!is.na(freq)]

p25 <- ggplot(fl, aes(reorder(feature, freq), freq, fill = arm)) +
  geom_col(position = "dodge") + coord_flip() +
  facet_wrap(~selector) +
  scale_y_continuous(labels = percent) +
  scale_fill_manual(values = c(control = "grey55", perturbed = "#C62828")) +
  labs(title = sprintf("Fig 25. How often each feature enters the top-%d subset", K_MAIN),
       subtitle = "Grey: control arm (fixed composition). Red: all composition-perturbed configurations.\nFeatures whose bars differ are the ones whose apparent importance depends on the attack mix.",
       x = NULL, y = "selection frequency", fill = NULL) +
  theme_ids()
save_fig(p25, "25_feature_frequency", w = 11, h = 9)
fwrite(freq_dt, file.path(REP_DIR, "07_feature_frequency.csv"))

saveRDS(list(rank_l1 = rank_l1, rank_rf = rank_rf, diag = diag,
             features = FEATS, k_values = K_VALUES, k_main = K_MAIN,
             stability = stab, jaccard = jac, seed = SEED,
             num_trees = NUM_TREES),
        file.path(OUT_DIR, "processed", "selection.rds"))

si <- c(
  sprintf("timestamp        : %s", format(Sys.time())),
  sprintf("R version        : %s", R.version.string),
  sprintf("glmnet / ranger  : %s / %s",
          as.character(packageVersion("glmnet")),
          as.character(packageVersion("ranger"))),
  sprintf("seed             : %d", SEED),
  sprintf("configurations   : %d", nrow(diag)),
  sprintf("selectors        : L1 logistic (glmnet), RF impurity (ranger, %d trees)", NUM_TREES),
  sprintf("k values         : %s", paste(K_VALUES, collapse = ", ")),
  sprintf("CONTROL Nogueira L1 (k=%d): %.4f", K_MAIN,
          ctrl[selector == "L1" & k == K_MAIN, nogueira]),
  sprintf("CONTROL Nogueira RF (k=%d): %.4f", K_MAIN,
          ctrl[selector == "RF" & k == K_MAIN, nogueira])
)
writeLines(si, file.path(REP_DIR, "07_run_record.txt"))
cat("\n", paste(si, collapse = "\n"), "\n")
cat("\nStage 7 complete. -> data/processed/selection.rds\n")
