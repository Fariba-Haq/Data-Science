

SEED    <- 13
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


sp <- readRDS(file.path(OUT_DIR, "processed", "split_safe.rds"))
y  <- sp$y8_train
n_train <- length(y)

before <- as.data.table(table(y)); setnames(before, c("class", "n"))
setorder(before, -n)
cat(sprintf("Stage 5: %s training rows, %d classes\n",
            format(n_train, big.mark = ","), nrow(before)))
print(before)
cat(sprintf("imbalance ratio (largest : smallest) = %.1f : 1\n",
            max(before$n) / min(before$n)))


resample_to <- function(y, k, seed) {
  set.seed(seed)
  unlist(lapply(sort(unique(y)), function(cl) {
    w <- which(y == cl)
    sample(w, k, replace = length(w) < k)
  }), use.names = FALSE)
}

k_under  <- min(before$n)
k_over   <- max(before$n)
k_hybrid <- as.integer(round(median(before$n)))

cat(sprintf("\ntarget per class:  under = %s | hybrid = %s | over = %s\n",
            comma(k_under), comma(k_hybrid), comma(k_over)))

idx <- list(
  undersample = resample_to(y, k_under,  SEED),
  hybrid      = resample_to(y, k_hybrid, SEED),
  oversample  = resample_to(y, k_over,   SEED)
)


cost <- rbindlist(lapply(names(idx), function(m) {
  ii <- idx[[m]]
  tb <- as.data.table(table(y[ii])); setnames(tb, c("class", "n_after"))
  tb <- merge(before, tb, by = "class")
  tb[, `:=`(method = m,
            duplication_factor = round(n_after / n, 2),
            unique_rows_used   = pmin(n, n_after))]
  tb[]
}))
cost[, method := factor(method, levels = c("undersample", "hybrid", "oversample"))]
setorder(cost, method, -n)
print(cost[, .(method, class, n_before = n, n_after,
               duplication_factor, unique_rows_used)], nrows = 30)

summ <- cost[, .(total_rows      = sum(n_after),
                 unique_rows     = sum(unique_rows_used),
                 rows_discarded  = sum(pmax(n - n_after, 0)),
                 max_duplication = max(duplication_factor)), by = method]
cat("\n--- cost of each strategy ---\n"); print(summ)
fwrite(cost, file.path(REP_DIR, "05_imbalance_cost.csv"))
fwrite(summ, file.path(REP_DIR, "05_imbalance_summary.csv"))



saveRDS(list(indices = idx,
             targets = c(undersample = k_under, hybrid = k_hybrid,
                         oversample = k_over),
             note = paste("Row indices into X_train / y8_train of",
                          "data/processed/split_safe.rds. Training partition only."),
             seed = SEED),
        file.path(OUT_DIR, "processed", "balanced_idx.rds"))
cat("\nsaved -> data/processed/balanced_idx.rds\n")


plot_dt <- rbind(
  cost[, .(method, class, n = n_before <- n, stage = "Before")],
  cost[, .(method, class, n = n_after,       stage = "After")]
)
plot_dt <- unique(plot_dt)
plot_dt[, stage := factor(stage, levels = c("Before", "After"))]

p18 <- ggplot(plot_dt, aes(reorder(class, n), n, fill = stage)) +
  geom_col(position = "dodge") + coord_flip() +
  facet_wrap(~method, nrow = 1) +
  scale_y_log10(labels = comma) +
  scale_fill_manual(values = c(Before = "grey55", After = "#2E7D32")) +
  labs(title = "Fig 18. Class balance before and after resampling",
       subtitle = sprintf("Applied to the training partition only - never before the split, never to the test set.\nLog scale. Imbalance before resampling is %.0f:1 (Benign : BruteForce).",
                          max(before$n) / min(before$n)),
       x = NULL, y = "training rows (log10)", fill = NULL) +
  theme_ids()
save_fig(p18, "18_imbalance", w = 12, h = 5.5)



si <- c(
  sprintf("timestamp        : %s", format(Sys.time())),
  sprintf("R version        : %s", R.version.string),
  sprintf("seed             : %d", SEED),
  sprintf("training rows    : %d", n_train),
  sprintf("classes          : %d", nrow(before)),
  sprintf("imbalance ratio  : %.1f : 1", max(before$n) / min(before$n)),
  sprintf("undersample к    : %d  (total %d rows)", k_under,  k_under  * nrow(before)),
  sprintf("hybrid target    : %d  (total %d rows)", k_hybrid, k_hybrid * nrow(before)),
  sprintf("oversample к     : %d  (total %d rows)", k_over,   k_over   * nrow(before)),
  sprintf("max duplication  : %.2fx (%s)",
          max(cost$duplication_factor),
          cost[which.max(duplication_factor), as.character(class)])
)
writeLines(si, file.path(REP_DIR, "05_run_record.txt"))
cat("\n", paste(si, collapse = "\n"), "\n")
cat("\nStage 5 complete. -> data/processed/balanced_idx.rds\n")
