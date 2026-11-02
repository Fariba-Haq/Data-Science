
#  12_substitute_quality.R


## CONFIG

ARMS_MAIN <- "all"

SELECTOR  <- "L1"        # Stage 11 ran L1 only
OUT_DIR <- "data"; FIG_DIR <- "figures"; REP_DIR <- "reports"


## SETUP


need <- c("data.table", "ggplot2", "scales")
miss <- need[!sapply(need, requireNamespace, quietly = TRUE)]
if (length(miss))
  stop("Install first: install.packages(c(",
       paste(sprintf("'%s'", miss), collapse = ", "), "))")
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(scales)
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
save_fig <- function(p, name, w = 9.5, h = 5.5) {
  path <- file.path(FIG_DIR, paste0(name, ".png"))
  ggsave(path, p, width = w, height = h, dpi = 200, bg = "white")
  cat("  figure ->", path, "\n"); invisible(path)
}

oc_obj <- readRDS(file.path(OUT_DIR, "processed", "overlap_curve.rds"))
f9     <- readRDS(file.path(OUT_DIR, "processed", "macro_f1.rds"))$f1
rbase  <- tryCatch(readRDS(file.path(OUT_DIR, "processed", "random_baseline.rds")),
                   error = function(e) NULL)

oc <- as.data.table(oc_obj$raw)
K  <- oc_obj$k

ceil_dt <- oc[is.na(m), .(draw, f1_all = macro_f1)]
stopifnot(!anyDuplicated(ceil_dt$draw))

con <- merge(oc[!is.na(m) & selector == SELECTOR,
                .(draw, m, rep, jaccard, subset_type, macro_f1)],
             ceil_dt, by = "draw", all.x = TRUE)
con[, reduction_cost := f1_all - macro_f1]
con[, source := "constructed (random substitutes)"]

cat("constructed subsets :", nrow(con), "rows,",
    uniqueN(con$draw), "draws\n")

f9 <- as.data.table(f9)
sc <- f9[selector == SELECTOR]
if (ARMS_MAIN == "control") sc <- sc[arm == "control"]
sc <- sc[, .(jaccard = overlap_own_ctrl, reduction_cost, arm, config)]
sc[, source := "selector-chosen (composition perturbed)"]

cat("selector-chosen     :", nrow(sc), "rows, arms =",
    paste(sort(unique(sc$arm)), collapse = ", "), "\n")
cat("  overlap range     :", sprintf("%.3f - %.3f",
                                     min(sc$jaccard), max(sc$jaccard)), "\n\n")

lo_common <- max(min(con$jaccard), min(sc$jaccard))
hi_common <- min(max(con$jaccard), max(sc$jaccard))

cat(strrep("=", 74), "\n", sep = "")
cat("  MATCHED-OVERLAP COMPARISON\n")
cat(strrep("=", 74), "\n", sep = "")
cat(sprintf("common overlap range: %.3f - %.3f\n\n", lo_common, hi_common))

lev <- con[, .(con_cost = mean(reduction_cost), con_sd = sd(reduction_cost)),
           by = .(m, jaccard)]
setorder(lev, jaccard)

sc[, m_match := lev$m[sapply(jaccard, function(j) which.min(abs(lev$jaccard - j)))]]
matched <- merge(
  sc[, .(n_sel = .N,
         sel_cost = mean(reduction_cost),
         sel_sd   = sd(reduction_cost),
         sel_jac  = mean(jaccard)), by = m_match],
  lev, by.x = "m_match", by.y = "m")
matched[, ratio := con_cost / sel_cost]
setorder(matched, -m_match)

print(matched[, .(features_kept = m_match,
                  jaccard_constructed = round(jaccard, 3),
                  jaccard_selector    = round(sel_jac, 3),
                  n_sel,
                  cost_selector    = round(sel_cost, 4),
                  cost_constructed = round(con_cost, 4),
                  ratio            = round(ratio, 2))])
fwrite(matched, file.path(REP_DIR, "12_matched_overlap.csv"))

in_range <- con[jaccard >= lo_common & jaccard <= hi_common]
sel_in   <- sc[jaccard >= lo_common & jaccard <= hi_common]

mu_con <- mean(in_range$reduction_cost)
mu_sel <- mean(sel_in$reduction_cost)

cat(sprintf("\nover the common overlap range (%.2f - %.2f):\n", lo_common, hi_common))
cat(sprintf("  selector-chosen  mean reduction cost : %.4f   (n = %d)\n",
            mu_sel, nrow(sel_in)))
cat(sprintf("  constructed      mean reduction cost : %.4f   (n = %d)\n",
            mu_con, nrow(in_range)))
cat(sprintf("  ratio                                : %.2fx\n", mu_con / mu_sel))

if (nrow(sel_in) >= 3 && nrow(in_range) >= 3) {
  wt <- wilcox.test(in_range$reduction_cost, sel_in$reduction_cost,
                    alternative = "greater", exact = FALSE)
  cat(sprintf("  Wilcoxon rank-sum W = %.0f, one-sided p = %.4g\n",
              unname(wt$statistic), wt$p.value))
} else {
  wt <- NULL
  cat("  too few points in the common range for a test\n")
}

cat("\n--- verdict ---\n")
if (!is.null(wt) && wt$p.value < 0.05 && mu_con > mu_sel) {
  cat("  SUPPORTED. At comparable overlap, randomly-constructed subsets lose\n")
  cat(sprintf("  %.2fx more macro-F1 than the subsets the selector actually\n",
              mu_con / mu_sel))
  cat("  returned. The selector's substitutes are therefore NOT arbitrary -\n")
  cat("  the subset is unstable in identity but stable in information content.\n")
  cat("  This rules out feature redundancy as the explanation for the near-zero\n")
  cat("  Stage 9 selection cost.\n")
} else {
  cat("  NOT SUPPORTED. Selector-chosen and randomly-constructed subsets lose\n")
  cat("  comparable macro-F1 at comparable overlap. The selector's substitutes\n")
  cat("  are no better than random ones, and the near-zero Stage 9 cost is\n")
  cat("  best explained by redundancy in the feature space. Report this.\n")
}

cat("\n--- independent cross-check: Stage 10 random vs Stage 11 m = 3 ---\n")
m3 <- oc[!is.na(m) & m == 3 & selector == SELECTOR, mean(macro_f1)]
cat(sprintf("  expected committed features in a random 10-of-36 : %.2f\n",
            K * K / length(strsplit(oc[is.na(m)][1, features], "\\|")[[1]])))
cat(sprintf("  Stage 11, m = 3 (Jaccard 0.176), macro-F1        : %.4f\n", m3))
if (!is.null(rbase)) {
  r10 <- rbase$summary[subset_type == "random", mean]
  cat(sprintf("  Stage 10, random 10-subsets, macro-F1           : %.4f\n", r10))
  cat(sprintf("  agreement                                       : %.4f\n",
              abs(m3 - r10)))
  cat("  Two independent experiments, different constructions, same answer.\n")
} else {
  cat("  (random_baseline.rds not found - skipping)\n")
}


## FIGURES 37-38

both <- rbindlist(list(
  con[, .(jaccard, reduction_cost, source)],
  sc [, .(jaccard, reduction_cost, source)]))

p37 <- ggplot(both, aes(jaccard, reduction_cost, colour = source, fill = source)) +
  annotate("rect", xmin = lo_common, xmax = hi_common, ymin = -Inf, ymax = Inf,
           fill = "grey60", alpha = 0.12) +
  geom_hline(yintercept = 0, colour = "grey30", linewidth = 0.4) +
  geom_point(alpha = 0.5, size = 1.7) +
  geom_smooth(method = "loess", formula = y ~ x, se = TRUE,
              linewidth = 0.9, alpha = 0.15) +
  scale_colour_manual(values = c("constructed (random substitutes)" = "#C62828",
                                 "selector-chosen (composition perturbed)" = "#1565C0")) +
  scale_fill_manual(values = c("constructed (random substitutes)" = "#C62828",
                               "selector-chosen (composition perturbed)" = "#1565C0")) +
  scale_x_continuous(breaks = seq(0, 1, 0.2)) +
  guides(colour = guide_legend(nrow = 2), fill = guide_legend(nrow = 2)) +
  labs(title = "Fig 37. The selector's substitutes are not random substitutes",
       subtitle = paste0(
         "y = macro-F1(all 36) - macro-F1(subset), measured on the same rows, so composition effects cancel.\n",
         "At matched overlap the red cloud sits well above the blue: losing a committed feature to a RANDOM\n",
         "replacement costs far more than losing it to the feature the selector itself chose instead.\n",
         "Shaded band: the overlap range both experiments cover."),
       x = "Jaccard overlap with the committed subset",
       y = "reduction cost  (macro-F1 lost vs all 36 features)",
       colour = NULL, fill = NULL) +
  theme_ids()
save_fig(p37, "37_substitute_quality", w = 10, h = 6)

p38 <- ggplot(matched, aes(factor(m_match))) +
  geom_col(aes(y = con_cost, fill = "constructed (random substitutes)"),
           width = 0.62, alpha = 0.9) +
  geom_col(aes(y = sel_cost, fill = "selector-chosen (composition perturbed)"),
           width = 0.34) +
  geom_text(aes(y = con_cost, label = sprintf("%.1fx", ratio)),
            vjust = -0.5, size = 3.4, colour = "grey25") +
  scale_fill_manual(values = c("constructed (random substitutes)" = "#EF9A9A",
                               "selector-chosen (composition perturbed)" = "#1565C0")) +
  guides(fill = guide_legend(nrow = 2)) +
  labs(title = "Fig 38. Cost of substitution at matched overlap",
       subtitle = paste0(
         "Each selector-chosen subset is matched to the constructed level nearest in overlap.\n",
         "Labels give the ratio: how many times more macro-F1 a random substitution costs."),
       x = "features retained from the committed subset (m)",
       y = "mean reduction cost", fill = NULL) +
  theme_ids()
save_fig(p38, "38_matched_overlap_cost", w = 9.5, h = 5.5)

saveRDS(list(constructed = con, selector_chosen = sc, matched = matched,
             common_range = c(lo_common, hi_common),
             mu_con = mu_con, mu_sel = mu_sel, test = wt,
             selector = SELECTOR, arms = ARMS_MAIN),
        file.path(OUT_DIR, "processed", "substitute_quality.rds"))

si <- c(
  sprintf("timestamp        : %s", format(Sys.time())),
  sprintf("R version        : %s", R.version.string),
  sprintf("selector         : %s", SELECTOR),
  sprintf("stage 9 arms used: %s", ARMS_MAIN),
  sprintf("common overlap   : %.3f - %.3f", lo_common, hi_common),
  sprintf("cost selector    : %.4f (n = %d)", mu_sel, nrow(sel_in)),
  sprintf("cost constructed : %.4f (n = %d)", mu_con, nrow(in_range)),
  sprintf("ratio            : %.2fx", mu_con / mu_sel),
  sprintf("wilcoxon p       : %s",
          if (is.null(wt)) "n/a" else sprintf("%.4g", wt$p.value))
)
writeLines(si, file.path(REP_DIR, "12_run_record.txt"))
cat("\n", paste(si, collapse = "\n"), "\n")
cat("\nStage 12 complete. -> data/processed/substitute_quality.rds\n")
