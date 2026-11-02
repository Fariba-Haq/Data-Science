
RAW_DIR <- "D:/0.1.Researches/DS/1.Dataset_author"
SEED    <- 13

OUT_DIR <- "data"
FIG_DIR <- "figures"
REP_DIR <- "reports"

CONSTANT_THRESH <- 0.999   
CORR_DUP_THRESH <- 0.999   
PLOT_SAMPLE     <- 200000 
VERIFY_ON_FULL  <- TRUE    


suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(scales)
})
for (d in c(FIG_DIR, REP_DIR)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

theme_ids <- function(base = 11) {
  theme_minimal(base_size = base) +
    theme(panel.grid.minor = element_blank(),
          plot.title    = element_text(face = "bold", size = base + 2),
          plot.subtitle = element_text(colour = "grey35", size = base - 1),
          strip.text    = element_text(face = "bold"),
          legend.position = "bottom")
}
save_fig <- function(p, name, w = 9, h = 5.5) {
  path <- file.path(FIG_DIR, paste0(name, ".png"))
  ggsave(path, p, width = w, height = h, dpi = 200, bg = "white")
  cat("  figure ->", path, "\n"); invisible(path)
}

inv     <- readRDS(file.path(OUT_DIR, "interim", "inventory.rds"))
ws      <- readRDS(file.path(OUT_DIR, "interim", "working_set.rds"))
all_csv <- inv$all_csv

feats <- setdiff(names(ws), c("label_fine", "label_8", "label_binary"))
Nw    <- nrow(ws)
cat(sprintf("Stage 3: %s rows x %d raw features\n",
            format(Nw, big.mark = ","), length(feats)))


col_stats <- rbindlist(lapply(feats, function(cn) {
  x   <- ws[[cn]]
  fin <- is.finite(x)
  xf  <- x[fin]
  u   <- unique(xf)
  data.table(
    column      = cn,
    n_unique    = length(u),
    sd          = if (length(xf) > 1) sd(xf) else 0,
    min         = if (length(xf)) min(xf) else NA_real_,
    max         = if (length(xf)) max(xf) else NA_real_,
    zero_share  = sum(xf == 0) / Nw,
    modal_share = if (!length(u)) NA_real_ else max(tabulate(match(xf, u))) / Nw,
    na_share    = sum(is.na(x)) / Nw,
    inf_share   = sum(is.infinite(x)) / Nw
  )
}))
col_stats[, zero_variance := (sd == 0) | (n_unique <= 1)]
col_stats[, near_constant := modal_share >= CONSTANT_THRESH]
setorder(col_stats, modal_share)

print(col_stats[, .(column, n_unique, sd = signif(sd, 4),
                    modal_share = round(modal_share, 4),
                    na_share = round(na_share, 6), zero_variance, near_constant)],
      nrows = 45)
cat("\nzero-variance:", sum(col_stats$zero_variance),
    "| near-constant at", CONSTANT_THRESH, ":", sum(col_stats$near_constant), "\n")
print(col_stats[na_share  > 0, .(column, n_na  = round(na_share  * Nw))])
print(col_stats[inf_share > 0, .(column, n_inf = round(inf_share * Nw))])
fwrite(col_stats, file.path(REP_DIR, "03_column_stats.csv"))


col_stats[, sig := paste(n_unique, signif(sd, 12), signif(min, 12), signif(max, 12))]
cand <- col_stats[, .N, by = sig][N > 1]

exact <- data.table(col_a = character(), col_b = character())
for (s in cand$sig) {
  g <- col_stats[sig == s, column]
  for (a in seq_along(g)) for (b in seq_along(g)) if (a < b)
    if (identical(ws[[g[a]]], ws[[g[b]]]))
      exact <- rbind(exact, data.table(col_a = g[a], col_b = g[b]))
}
cat("\nexact duplicate pairs:\n"); print(exact)


set.seed(SEED)
samp <- ws[sample.int(.N, min(PLOT_SAMPLE, .N))]
M <- as.matrix(samp[, ..feats]); M[!is.finite(M)] <- NA
cm <- suppressWarnings(cor(M, use = "pairwise.complete.obs"))
cm[!is.finite(cm)] <- 0
rm(M); gc(verbose = FALSE)

cdt <- as.data.table(as.table(cm)); setnames(cdt, c("f1", "f2", "r"))
cdt <- cdt[as.character(f1) < as.character(f2)]
cat("\n|r| >=", CORR_DUP_THRESH, ":\n");  print(cdt[abs(r) >= CORR_DUP_THRESH][order(-abs(r))])
cat("\n0.95 <= |r| <", CORR_DUP_THRESH, ":\n")
print(cdt[abs(r) >= 0.95 & abs(r) < CORR_DUP_THRESH][order(-abs(r))])
fwrite(cdt[abs(r) >= 0.95], file.path(REP_DIR, "03_correlated_pairs.csv"))
saveRDS(cm, file.path(OUT_DIR, "interim", "cor_matrix.rds"))


BORDERLINE <- col_stats[modal_share >= 0.98][order(-modal_share), column]
lab <- ws$label_8
res <- rbindlist(lapply(BORDERLINE, function(cn) {
  d <- data.table(label_8 = lab, nz = ws[[cn]] > 0)
  d[, .(n = .N, n_nz = sum(nz, na.rm = TRUE)), by = label_8][, column := cn][]
}))
res[, nz_rate     := n_nz / n]
res[, share_of_nz := n_nz / sum(n_nz), by = column]

cat("\n% of each category's rows where the feature is non-zero:\n")
w1 <- dcast(res, column ~ label_8, value.var = "nz_rate")
num <- setdiff(names(w1), "column")
w1[, (num) := lapply(.SD, function(x) round(100 * x, 3)), .SDcols = num]
print(w1)
fwrite(res, file.path(REP_DIR, "03_borderline_by_category.csv"))


DROP_DUPLICATE <- c("Tot size", "LLC", "IPv")
DROP_CONSTANT  <- character(0)          # set to c("Telnet","SMTP") to drop them
DROP_COLS      <- c(DROP_DUPLICATE, DROP_CONSTANT)
KEEP_FEATS     <- setdiff(feats, DROP_COLS)

cat(sprintf("\nraw %d - dropped %d = FINAL %d features\n",
            length(feats), length(DROP_COLS), length(KEEP_FEATS)))

sens <- data.table(threshold = c(0.95, 0.99, 0.995, 0.997, 0.999, 0.9995))
sens[, n_near_constant := sapply(threshold, function(t) sum(col_stats$modal_share >= t))]
sens[, features_left   := length(feats) - length(DROP_DUPLICATE) - n_near_constant]
print(sens)

if (VERIFY_ON_FULL) {
  chk_cols <- c("AVG", "Tot size", "IPv", "LLC", "ARP")
  d1 <- d2 <- d3 <- 0
  for (i in seq_along(all_csv)) {
    d <- fread(all_csv[i], select = chk_cols, showProgress = FALSE)
    d1 <- max(d1, max(abs(d[["AVG"]] - d[["Tot size"]])))
    d2 <- max(d2, max(abs(d[["IPv"]] - d[["LLC"]])))
    d3 <- max(d3, max(abs(d[["ARP"]] - (1 - d[["IPv"]]))))
    rm(d)
    if (i %% 100 == 0 || i == length(all_csv)) cat(sprintf("  verify [%d/%d]\n", i, length(all_csv)))
  }
  cat(sprintf("full-data check: |AVG-Tot size|=%g  |IPv-LLC|=%g  |ARP-(1-IPv)|=%g\n",
              d1, d2, d3))
}

rule <- c(
  "CLEANING RULE (as applied)",
  "",
  sprintf("Raw feature columns in this release : %d", length(feats)),
  "Removed - exact duplicates / perfect collinearity:",
  "  Tot size == AVG      (bit-identical, all 46,776,697 rows)",
  "  LLC      == IPv      (bit-identical, all 46,776,697 rows)",
  "  IPv      == 1 - ARP  (r = -1.000, max deviation 2.5e-16)",
  sprintf("Retained                            : %d", length(KEEP_FEATS)),
  "",
  "Near-constant columns were NOT removed. Sparsity in this dataset is",
  "class-conditional (SSH: 32.1% non-zero in BruteForce vs <=2.6% elsewhere;",
  "Telnet/SMTP: non-zero only when Recon/DoS/DDoS are present). Removing such",
  "columns would confound the composition-sensitivity experiment.",
  "",
  "NOTE: 36 here is 39 raw minus cleaning. Keskin & Oz and Aragao also report",
  "39, but theirs is 47 minus cleaning. Different route, similar number - say",
  "so explicitly in Methods.",
  "",
  "Data-quality notes:",
  "  Std, Variance : 100 NA rows   (4.66e-05)",
  "  Rate          : 104 Inf rows  (4.84e-05)",
  "  3 source files have a truncated final line, discarded by fread.",
  "  All handled by train-fold median imputation in Stage 4."
)
writeLines(rule, file.path(REP_DIR, "03_cleaning_rule.txt"))
saveRDS(KEEP_FEATS, file.path(OUT_DIR, "interim", "keep_features.rds"))
cat("\n", paste(rule, collapse = "\n"), "\n")



## Fig 7 
cs <- copy(col_stats)[, status := ifelse(column %in% DROP_COLS,
                                         "dropped (redundant)", "retained")]
p7 <- ggplot(cs, aes(reorder(column, modal_share), modal_share, fill = status)) +
  geom_col() + coord_flip() +
  geom_hline(yintercept = CONSTANT_THRESH, linetype = "dashed", colour = "#C62828") +
  annotate("text", x = 3, y = 0.985, label = "0.999 threshold  ",
           hjust = 1, size = 3, colour = "#C62828") +
  scale_y_continuous(labels = percent, limits = c(0, 1.02)) +
  scale_fill_manual(values = c(retained = "grey60",
                               `dropped (redundant)` = "#C62828")) +
  labs(title = "Fig 7. Modal-value share per feature",
       subtitle = "No column has zero variance. The three dropped columns were removed as exact duplicates,\nnot for near-constancy - see Fig 8 for why the constancy threshold was not applied.",
       x = NULL, y = "share of rows holding the modal value", fill = NULL) +
  theme_ids()
save_fig(p7, "07_feature_constancy", h = 9)

## Fig 8
p8 <- ggplot(res, aes(label_8, reorder(column, nz_rate), fill = 100 * nz_rate)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%.2f", 100 * nz_rate)), size = 3, colour = "grey20") +
  scale_fill_gradient(low = "white", high = "#C62828", name = "% non-zero") +
  labs(title = "Fig 8. Sparse features are class-conditional, not uninformative",
       subtitle = "SSH fires in 32.1% of BruteForce rows against <=2.6% elsewhere. Removing near-constant\ncolumns would delete signal whose visibility depends on attack composition.",
       x = NULL, y = NULL) +
  theme_ids() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1), panel.grid = element_blank())
save_fig(p8, "08_class_conditional_sparsity", w = 10, h = 5.5)

## Fig 9 
cmk  <- cm[KEEP_FEATS, KEEP_FEATS]
ord  <- hclust(as.dist(1 - abs(cmk)))$order     # distance = 1 - |r|
lv   <- KEEP_FEATS[ord]
cdt2 <- as.data.table(as.table(cmk)); setnames(cdt2, c("f1", "f2", "r"))
cdt2[, f1 := factor(as.character(f1), levels = lv)]
cdt2[, f2 := factor(as.character(f2), levels = lv)]

p9 <- ggplot(cdt2, aes(f1, f2, fill = r)) + geom_tile() +
  scale_fill_gradient2(low = "#1565C0", mid = "white", high = "#C62828",
                       midpoint = 0, limits = c(-1, 1)) +
  labs(title = sprintf("Fig 9. Correlation among the %d retained features", length(KEEP_FEATS)),
       subtitle = "Clustered ordering. Redundant blocks motivate selection; exact duplicates were already removed.",
       x = NULL, y = NULL, fill = "Pearson r") +
  theme_ids() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
        axis.text.y = element_text(size = 7), panel.grid = element_blank())
save_fig(p9, "09_correlation_heatmap", w = 10, h = 9)

## Fig 10 
skew1 <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 3) return(0)
  s <- sd(x); if (s == 0) return(0)
  mean((x - mean(x))^3) / s^3
}
sk <- data.table(column = KEEP_FEATS,
                 skew = sapply(KEEP_FEATS, function(c) skew1(ws[[c]])))
sk <- sk[order(-abs(skew))]            
print(sk[1:12]); fwrite(sk, file.path(REP_DIR, "03_skewness.csv"))

p10 <- ggplot(sk, aes(reorder(column, abs(skew)), skew, fill = abs(skew) > 5)) +
  geom_col() + coord_flip() +
  scale_fill_manual(values = c(`FALSE` = "grey60", `TRUE` = "#C62828"),
                    labels = c("moderate", "heavy-tailed"), name = NULL) +
  labs(title = "Fig 10. Skewness per retained feature",
       subtitle = "Heavy tails are why scaling uses median/MAD rather than mean/SD, with clipping at +/-3",
       x = NULL, y = "skewness") +
  theme_ids()
save_fig(p10, "10_skewness", h = 9)

## Fig 11 — distributions.
CONT_FEATS <- col_stats[n_unique >= 1000 & column %in% KEEP_FEATS][
                order(-n_unique)][1:min(12, .N), column]
cat("\ncontinuous features plotted in Fig 11:\n"); print(CONT_FEATS)

set.seed(SEED)
sub  <- ws[sample.int(.N, min(PLOT_SAMPLE, .N))]
long <- melt(sub[, c(CONT_FEATS, "label_binary"), with = FALSE],
             id.vars = "label_binary", variable.name = "feature", value.name = "value")
long <- long[is.finite(value)]
long[, cls   := ifelse(label_binary == 1, "Attack", "Benign")]
long[, vlog  := sign(value) * log10(1 + abs(value))]   # signed log1p

p11 <- ggplot(long, aes(vlog, fill = cls)) +
  geom_density(alpha = 0.5, colour = NA) +
  facet_wrap(~feature, scales = "free", ncol = 4) +
  scale_fill_manual(values = c(Benign = "#2E7D32", Attack = "#C62828")) +
  labs(title = "Fig 11. Distributions of the continuous features, before scaling",
       subtitle = "Signed log10(1+|x|) axis. Benign and attack separate visibly in several features,\nand the long right tails are what median/MAD scaling with clipping is there to handle.",
       x = "signed log10(1 + |value|)", y = "density", fill = NULL) +
  theme_ids(9)
save_fig(p11, "11_distributions_continuous", w = 12, h = 8)

cat("\nStage 3 complete. -> data/interim/keep_features.rds\n")
