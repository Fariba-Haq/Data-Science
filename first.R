##step-01

install.packages(c("data.table", "ggplot2", "scales"))

library(data.table)
library(ggplot2)
library(scales)

for (d in c("data", "data/interim", "data/processed", "figures", "reports")) {
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
}

cat("project dir:", getwd(), "\n")
list.dirs(".", recursive = FALSE)

RAW_DIR <- "D:\\0.1.Researches\\DS\\1.Dataset_author"   

if (!dir.exists(RAW_DIR)) stop("RAW_DIR not found. Check the path.")
cat("OK:", normalizePath(RAW_DIR), "\n")
head(list.files(RAW_DIR), 40)

folders <- list.dirs(RAW_DIR, recursive = FALSE, full.names = TRUE)

inv <- data.table(
  folder = basename(folders),
  n_csv  = sapply(folders, function(d)
    length(list.files(d, pattern = "\\.csv$", ignore.case = TRUE))),
  MB     = sapply(folders, function(d) {
    f <- list.files(d, pattern = "\\.csv$", full.names = TRUE,
                    ignore.case = TRUE)
    if (!length(f)) 0 else round(sum(file.info(f)$size) / 1e6, 1)
  })
)
setorder(inv, -MB)

print(inv, nrows = 40)

cat("\nfolders      :", nrow(inv),
    "\ntotal CSVs   :", sum(inv$n_csv),
    "\ntotal size   :", round(sum(inv$MB) / 1024, 2), "GB",
    "\nempty folders:", sum(inv$n_csv == 0), "\n")

fwrite(inv, "reports/00_folder_inventory.csv")

inv[, share_pct := round(100 * MB / sum(MB), 3)]
print(inv[, .(folder, n_csv, MB, share_pct)], nrows = 40)

cat("\n--- benign vs attack (by size, rough proxy for rows) ---\n")
print(inv[grepl("Benign", folder, ignore.case = TRUE)])

#6
all_csv <- list.files(RAW_DIR, pattern = "\\.csv$", full.names = TRUE,
                      recursive = TRUE, ignore.case = TRUE)
cat("files found:", length(all_csv), "\n")

hdrs <- lapply(all_csv, function(f) names(fread(f, nrows = 0)))
sig  <- sapply(hdrs, paste, collapse = "|")
u    <- unique(sig)

cat("distinct header signatures:", length(u), "\n\n")

if (length(u) == 1) {
  cat("GOOD - all 309 files share one schema.\n")
} else {
  cat("!! schemas differ. Files per variant:\n")
  print(table(match(sig, u)))
  for (v in seq_along(u)) {
    cat("\n--- variant", v, "---\n")
    print(basename(all_csv[sig == u[v]])[1:3])
    print(setdiff(strsplit(u[v], "|", fixed = TRUE)[[1]],
                  strsplit(u[1], "|", fixed = TRUE)[[1]]))
  }
}

ref <- hdrs[[1]]
cat("\ncolumns:", length(ref), "\n")
print(data.frame(n = seq_along(ref), column = ref))

cat("\nmagnitude-like column:",
    grep("magnit", ref, ignore.case = TRUE, value = TRUE), "\n")
cat("label-like column    :",
    grep("^label|^class$|^attack$", ref, ignore.case = TRUE, value = TRUE), "\n")

#7
f1 <- all_csv[1]
cat("file:", basename(f1), "\n\n")

peek <- fread(f1, nrows = 5)
print(str(peek))
cat("\n--- first 3 rows, first 12 cols ---\n")
print(peek[1:3, 1:12])

cat("\ncolumn classes:\n")
print(table(sapply(peek, function(x) class(x)[1])))

#8
KAGGLE_DIR <- "D:/0.1.Researches/DS/3.Dataset_kaggle"

cat("exists:", dir.exists(KAGGLE_DIR), "\n\n")
print(list.files(KAGGLE_DIR, recursive = TRUE)[1:20])

kcsv <- list.files(KAGGLE_DIR, pattern = "\\.csv$", full.names = TRUE,
                   recursive = TRUE, ignore.case = TRUE)
cat("\nkaggle csv files:", length(kcsv), "\n")

if (length(kcsv)) {
  kh <- names(fread(kcsv[1], nrows = 0))
  cat("kaggle columns:", length(kh), "\n")
  print(data.frame(n = seq_along(kh), column = kh))
  
  cat("\nin kaggle, NOT in author's 39:\n"); print(setdiff(kh, ref))
  cat("\nin author's 39, NOT in kaggle:\n"); print(setdiff(ref, kh))
}

#9
for (f in kcsv) {
  lab <- fread(f, select = "label", showProgress = FALSE)
  cat("\n=====", basename(f), "=====\n")
  cat("rows:", format(nrow(lab), big.mark = ","),
      " | size:", round(file.info(f)$size / 1e6, 1), "MB\n")
  cat("distinct labels:", uniqueN(lab$label), "\n")
  tb <- sort(table(lab$label), decreasing = TRUE)
  print(head(tb, 8))
  cat("benign share:",
      sprintf("%.3f%%", 100 * sum(grepl("Benign", lab$label, ignore.case = TRUE)) / nrow(lab)),
      "\n")
  rm(lab); gc(verbose = FALSE)
}

#10
t0 <- Sys.time()
res <- vector("list", length(all_csv))

for (i in seq_along(all_csv)) {
  n <- tryCatch(nrow(fread(all_csv[i], select = 1L, showProgress = FALSE)),
                error = function(e) NA_integer_)
  res[[i]] <- data.table(
    label = basename(dirname(all_csv[i])),
    file  = basename(all_csv[i]),
    rows  = n
  )
  if (i %% 25 == 0 || i == length(all_csv))
    cat(sprintf("[%d/%d] %s elapsed\n", i, length(all_csv),
                format(round(difftime(Sys.time(), t0)))))
}

rowcounts <- rbindlist(res)
fwrite(rowcounts, "reports/01_rows_per_file.csv")

cat("\nfailed files:", sum(is.na(rowcounts$rows)), "\n")
cat("total rows  :", format(sum(rowcounts$rows, na.rm = TRUE), big.mark = ","), "\n")

#11

lab_dt <- rowcounts[, .(n_files = .N, rows = sum(rows, na.rm = TRUE)), by = label]
lab_dt[, share_pct := round(100 * rows / sum(rows), 4)]
setorder(lab_dt, -rows)
print(lab_dt, nrows = 40)

N <- sum(lab_dt$rows)
benign <- lab_dt[grepl("Benign", label, ignore.case = TRUE), sum(rows)]

cat("\n--- HEADLINE NUMBERS (these go in the paper) ---\n")
cat("total rows   :", format(N, big.mark = ","), "\n")
cat("labels       :", nrow(lab_dt), "\n")
cat("benign rows  :", format(benign, big.mark = ","), "\n")
cat("benign share :", sprintf("%.4f%%", 100 * benign / N), "\n")
cat("smallest class:", lab_dt[which.min(rows), label],
    "=", format(min(lab_dt$rows), big.mark = ","), "rows\n")

fwrite(lab_dt, "reports/01_label_counts.csv")

#12
CAT_MAP <- c(
  # DDoS (12)
  "DDoS-ICMP_Flood"="DDoS", "DDoS-UDP_Flood"="DDoS", "DDoS-TCP_Flood"="DDoS",
  "DDoS-PSHACK_FLOOD"="DDoS", "DDoS-SYN_Flood"="DDoS", "DDoS-RSTFINFLOOD"="DDoS",
  "DDoS-SynonymousIP_Flood"="DDoS", "DDoS-ICMP_Fragmentation"="DDoS",
  "DDoS-UDP_Fragmentation"="DDoS", "DDoS-ACK_Fragmentation"="DDoS",
  "DDoS-HTTP_Flood"="DDoS", "DDoS-SlowLoris"="DDoS",
  # DoS (4)
  "DoS-UDP_Flood"="DoS", "DoS-TCP_Flood"="DoS",
  "DoS-SYN_Flood"="DoS", "DoS-HTTP_Flood"="DoS",
  # Mirai (3)
  "Mirai-greeth_flood"="Mirai", "Mirai-greip_flood"="Mirai",
  "Mirai-udpplain"="Mirai",
  # Recon (5)
  "VulnerabilityScan"="Recon", "Recon-HostDiscovery"="Recon",
  "Recon-OSScan"="Recon", "Recon-PortScan"="Recon", "Recon-PingSweep"="Recon",
  # Spoofing (2)
  "MITM-ArpSpoofing"="Spoofing", "DNS_Spoofing"="Spoofing",
  # Web (6)
  "BrowserHijacking"="Web", "CommandInjection"="Web", "SqlInjection"="Web",
  "XSS"="Web", "Backdoor_Malware"="Web", "Uploading_Attack"="Web",
  # BruteForce (1)
  "DictionaryBruteForce"="BruteForce",
  # Benign (1)
  "Benign_Final"="Benign"
)

# hard verification: every folder maps, every map entry exists on disk
stopifnot(length(CAT_MAP) == 34)
cat("folders with no map entry:",
    paste(setdiff(lab_dt$label, names(CAT_MAP)), collapse = ", "), "\n")
cat("map entries with no folder:",
    paste(setdiff(names(CAT_MAP), lab_dt$label), collapse = ", "), "\n")

BENIGN_FOLDER <- "Benign_Final"
ATTACK_CATS   <- c("DDoS","DoS","Mirai","Recon","Spoofing","Web","BruteForce")

#13
lab_dt[, category := CAT_MAP[label]]

cat_dt <- lab_dt[, .(n_labels = .N, n_files = sum(n_files), rows = sum(rows)),
                 by = category]
cat_dt[, share_pct := round(100 * rows / sum(rows), 4)]
setorder(cat_dt, -rows)
print(cat_dt)

fwrite(lab_dt, "reports/01_label_counts.csv")
fwrite(cat_dt, "reports/01_category_counts.csv")

# --- can Stage 6 build its configurations? ---
CONFIG_ROWS <- 100000; BENIGN_FRACTION <- 0.5; TEST_FRACTION <- 0.30
n_attack_needed <- CONFIG_ROWS * (1 - BENIGN_FRACTION)

feas <- cat_dt[category != "Benign", .(category, rows)]
feas[, train_rows := floor(rows * (1 - TEST_FRACTION))]
feas[, max_share_of_attack_half := round(pmin(1, train_rows / n_attack_needed), 3)]
setorder(feas, max_share_of_attack_half)
print(feas)

cat("\nbenign train rows available:",
    format(floor(benign * (1 - TEST_FRACTION)), big.mark = ","),
    " | needed per config:", format(CONFIG_ROWS * BENIGN_FRACTION, big.mark = ","), "\n")

#14
SEED            <- 13
CATEGORY_CAP    <- 200000
MIN_PER_LABEL   <- 3000
BENIGN_TAKE_ALL <- TRUE

# ---- label-level quotas ----
lab_dt[, cat_rows   := sum(rows), by = category]
lab_dt[, cat_target := pmin(cat_rows, CATEGORY_CAP)]
if (BENIGN_TAKE_ALL) lab_dt[category == "Benign", cat_target := cat_rows]

lab_dt[, quota := pmin(rows,
                       pmax(round(cat_target * rows / cat_rows),
                            pmin(MIN_PER_LABEL, rows)))]

print(lab_dt[order(category, -rows),
             .(label, category, rows, quota,
               pct_kept = round(100 * quota / rows, 2))], nrows = 40)

cat("\n--- working set by category ---\n")
plan_cat <- lab_dt[, .(rows_full = sum(rows), rows_kept = sum(quota)), by = category]
plan_cat[, pct_of_ws := round(100 * rows_kept / sum(rows_kept), 2)]
setorder(plan_cat, -rows_kept)
print(plan_cat)

WS_ROWS <- sum(lab_dt$quota)
cat("\nworking set rows :", format(WS_ROWS, big.mark = ","), "\n")
cat("est. memory      :", round(WS_ROWS * 39 * 8 / 1e9, 2), "GB (numeric cols only)\n")

# ---- file-level allocation, proportional to each file's true size ----
plan <- merge(rowcounts, lab_dt[, .(label, label_rows = rows, quota)], by = "label")
setorder(plan, label, -rows)
plan[, take := round(quota * rows / label_rows)]
plan[, take := { t <- take; t[1] <- t[1] + (quota[1] - sum(t)); t }, by = label]
plan[, take := pmax(pmin(take, rows), 0L)]

chk <- plan[, .(allocated = sum(take), quota = quota[1]), by = label][allocated != quota]
cat("\nlabels whose file allocation != quota:", nrow(chk), "\n")
if (nrow(chk)) print(chk)

cat("total rows to read:", format(sum(plan$take), big.mark = ","), "\n")
fwrite(plan, "reports/02_sampling_plan.csv")

#15
set.seed(SEED)
setorder(plan, label, -rows, file)      # fixed order => reproducible sampling

t0 <- Sys.time()
parts <- vector("list", nrow(plan))

for (i in seq_len(nrow(plan))) {
  if (plan$take[i] <= 0) next
  f <- file.path(RAW_DIR, plan$label[i], plan$file[i])
  
  dt <- fread(f, showProgress = FALSE)
  stopifnot(identical(names(dt), ref))          # schema guard
  
  if (plan$take[i] < nrow(dt))
    dt <- dt[sample.int(nrow(dt), plan$take[i])]
  
  dt[, label_fine := plan$label[i]]
  parts[[i]] <- dt
  rm(dt)
  
  if (i %% 25 == 0 || i == nrow(plan)) {
    gc(verbose = FALSE)
    cat(sprintf("[%d/%d] %s\n", i, nrow(plan),
                format(round(difftime(Sys.time(), t0)))))
  }
}

ws <- rbindlist(parts, use.names = TRUE)
rm(parts); gc(verbose = FALSE)

ws[, label_8      := CAT_MAP[label_fine]]
ws[, label_binary := as.integer(label_fine != BENIGN_FOLDER)]
ws <- ws[sample.int(.N)]                        # shuffle once, here

saveRDS(ws, "data/interim/working_set.rds")

cat("\nrows :", format(nrow(ws), big.mark = ","),
    "\ncols :", ncol(ws),
    "\nmemory:", round(as.numeric(object.size(ws)) / 1e9, 2), "GB\n")
cat("elapsed:", format(round(difftime(Sys.time(), t0))), "\n")

#16
# every label present at exactly its quota?
got <- ws[, .N, by = label_fine]
vfy <- merge(lab_dt[, .(label_fine = label, quota)], got, by = "label_fine")
vfy[, delta := N - quota]
cat("labels off quota:", sum(vfy$delta != 0), "\n")
if (any(vfy$delta != 0)) print(vfy[delta != 0])

cat("\nNA in label_8:", sum(is.na(ws$label_8)), "\n")
print(ws[, .N, by = label_8][order(-N)])
cat("\nbenign share of working set:",
    sprintf("%.2f%%", 100 * mean(ws$label_binary == 0)), "\n")

# reproducibility record — the brief asks for this per run
si <- c(
  sprintf("timestamp      : %s", format(Sys.time())),
  sprintf("R version      : %s", R.version.string),
  sprintf("platform       : %s", R.version$platform),
  sprintf("data.table     : %s", as.character(packageVersion("data.table"))),
  sprintf("seed           : %d", SEED),
  sprintf("CATEGORY_CAP   : %d", CATEGORY_CAP),
  sprintf("MIN_PER_LABEL  : %d", MIN_PER_LABEL),
  sprintf("BENIGN_TAKE_ALL: %s", BENIGN_TAKE_ALL),
  sprintf("source rows    : %d", sum(lab_dt$rows)),
  sprintf("working set    : %d", nrow(ws))
)
writeLines(si, "reports/02_run_record.txt")
cat("\n", paste(si, collapse = "\n"), "\n")

#17
ws <- readRDS("data/interim/working_set.rds")   # skip if already in session
feats <- setdiff(names(ws), c("label_fine", "label_8", "label_binary"))
Nw <- nrow(ws)

col_stats <- rbindlist(lapply(feats, function(cn) {
  x <- ws[[cn]]
  fin <- is.finite(x)
  xf  <- x[fin]
  u   <- unique(xf)
  modal <- if (length(u) == 0) NA_real_ else max(tabulate(match(xf, u))) / Nw
  data.table(
    column      = cn,
    n_unique    = length(u),
    sd          = if (length(xf) > 1) sd(xf) else 0,
    min         = if (length(xf)) min(xf) else NA_real_,
    max         = if (length(xf)) max(xf) else NA_real_,
    zero_share  = sum(xf == 0) / Nw,
    modal_share = modal,
    na_share    = sum(is.na(x)) / Nw,
    inf_share   = sum(is.infinite(x)) / Nw
  )
}))

CONSTANT_THRESH <- 0.999
col_stats[, zero_variance := (sd == 0) | (n_unique <= 1)]
col_stats[, near_constant := modal_share >= CONSTANT_THRESH]

setorder(col_stats, modal_share)
print(col_stats[, .(column, n_unique, sd = signif(sd, 4),
                    zero_share = round(zero_share, 4),
                    modal_share = round(modal_share, 4),
                    na_share = round(na_share, 6),
                    zero_variance, near_constant)], nrows = 45)

cat("\nzero-variance columns:", sum(col_stats$zero_variance), "\n")
cat("near-constant columns:", sum(col_stats$near_constant), "\n")
cat("columns with any NA  :", sum(col_stats$na_share > 0), "\n")
cat("columns with any Inf :", sum(col_stats$inf_share > 0), "\n")

fwrite(col_stats, "reports/03_column_stats.csv")

#18
# --- which columns carry NA / Inf ---
cat("columns with NA:\n");  print(col_stats[na_share  > 0, .(column, na_share,  n_bad = round(na_share  * Nw))])
cat("\ncolumns with Inf:\n"); print(col_stats[inf_share > 0, .(column, inf_share, n_bad = round(inf_share * Nw))])

# --- exact duplicates: group by statistical signature, then confirm ---
col_stats[, sig := paste(n_unique, signif(sd, 12), signif(min, 12), signif(max, 12))]
cand <- col_stats[, .N, by = sig][N > 1]

exact <- data.table(col_a = character(), col_b = character())
for (s in cand$sig) {
  g <- col_stats[sig == s, column]
  for (a in seq_along(g)) for (b in seq_along(g)) if (a < b)
    if (identical(ws[[g[a]]], ws[[g[b]]]))
      exact <- rbind(exact, data.table(col_a = g[a], col_b = g[b]))
}
cat("\nEXACT duplicate pairs:\n"); print(exact)

# --- near-duplicates by correlation, on a 200k sample ---
set.seed(SEED)
samp <- ws[sample.int(.N, 200000), ..feats]
M <- as.matrix(samp); M[!is.finite(M)] <- NA
cm <- suppressWarnings(cor(M, use = "pairwise.complete.obs"))
cm[!is.finite(cm)] <- 0

cdt <- as.data.table(as.table(cm)); setnames(cdt, c("f1", "f2", "r"))
cdt <- cdt[as.character(f1) < as.character(f2)]
cat("\npairs with |r| >= 0.999:\n")
print(cdt[abs(r) >= 0.999][order(-abs(r))])
cat("\npairs with 0.95 <= |r| < 0.999:\n")
print(cdt[abs(r) >= 0.95 & abs(r) < 0.999][order(-abs(r))])

rm(M, samp); gc(verbose = FALSE)
fwrite(cdt[abs(r) >= 0.95], "reports/03_correlated_pairs.csv")
saveRDS(cm, "data/interim/cor_matrix.rds")   # reused for the Fig 9 heatmap

#19

# confirm the ARP/IPv relationship exactly, on all rows
cat("ARP == 1 - IPv exactly:", identical(ws$ARP, 1 - ws$IPv), "\n")
cat("max |ARP - (1 - IPv)| :", max(abs(ws$ARP - (1 - ws$IPv))), "\n\n")

BORDERLINE <- col_stats[modal_share >= 0.98][order(-modal_share), column]
cat("testing:", paste(BORDERLINE, collapse = ", "), "\n\n")

lab <- ws$label_8
res <- rbindlist(lapply(BORDERLINE, function(cn) {
  d <- data.table(label_8 = lab, nz = ws[[cn]] > 0)
  out <- d[, .(n = .N, n_nz = sum(nz, na.rm = TRUE)), by = label_8]
  out[, column := cn][]
}))
res[, nz_rate      := n_nz / n]
res[, share_of_nz  := n_nz / sum(n_nz), by = column]

cat("--- % of each category's rows where the feature is NON-ZERO ---\n")
w1 <- dcast(res, column ~ label_8, value.var = "nz_rate")
num <- setdiff(names(w1), "column")
w1[, (num) := lapply(.SD, function(x) round(100 * x, 3)), .SDcols = num]
print(w1)

cat("\n--- where each feature's non-zero mass LIVES (% of its non-zero rows) ---\n")
w2 <- dcast(res, column ~ label_8, value.var = "share_of_nz")
num2 <- setdiff(names(w2), "column")
w2[, (num2) := lapply(.SD, function(x) round(100 * x, 2)), .SDcols = num2]
print(w2)

cat("\n--- concentration: top category and its share of non-zero mass ---\n")
conc <- res[, .SD[which.max(share_of_nz)], by = column][
  , .(column, top_category = label_8,
      pct_of_nonzero = round(100 * share_of_nz, 2),
      n_nonzero = n_nz)]
setorder(conc, -pct_of_nonzero)
print(conc)

fwrite(res, "reports/03_borderline_by_category.csv")

#20

DROP_DUPLICATE <- c("Tot size", "LLC", "IPv")   # AVG kept; ARP kept
DROP_CONSTANT  <- character(0)                  # <- set to c("Telnet","SMTP") to drop them instead

DROP_COLS  <- c(DROP_DUPLICATE, DROP_CONSTANT)
KEEP_FEATS <- setdiff(feats, DROP_COLS)

cat("raw features   :", length(feats), "\n")
cat("dropped        :", length(DROP_COLS), "->", paste(DROP_COLS, collapse = ", "), "\n")
cat("FINAL FEATURES :", length(KEEP_FEATS), "\n\n")

# threshold sensitivity — report this, don't hide it
sens <- data.table(threshold = c(0.95, 0.99, 0.995, 0.997, 0.999, 0.9995))
sens[, n_near_constant := sapply(threshold, function(t) sum(col_stats$modal_share >= t))]
sens[, features_left := length(feats) - length(DROP_DUPLICATE) - n_near_constant]
print(sens)

rule <- c(
  "CLEANING RULE (as applied)",
  "",
  sprintf("Raw feature columns in this release        : %d", length(feats)),
  "Removed, exact duplicates / perfect collinearity:",
  "  Tot size  == AVG            (bit-identical)",
  "  LLC       == IPv            (bit-identical)",
  "  IPv       == 1 - ARP        (r = -1.000, max dev 2.4e-16)",
  sprintf("Retained                                    : %d", length(KEEP_FEATS)),
  "",
  "Near-constant columns were NOT removed. Justification:",
  "  Sparsity in this dataset is class-conditional (e.g. SSH is non-zero in",
  "  32.1% of BruteForce rows vs <=2.6% elsewhere; Telnet is non-zero only",
  "  when Recon/DoS/DDoS are present). Removing such columns during",
  "  preprocessing would confound the composition-sensitivity experiment.",
  "",
  "Data-quality notes:",
  "  Std, Variance : 100 NA rows   (4.66e-05)",
  "  Rate          : 104 Inf rows  (4.84e-05)",
  "  3 source files have a truncated final line, discarded by fread (3 rows).",
  "  Handled by train-fold median imputation in Stage 4."
)
writeLines(rule, "reports/03_cleaning_rule.txt")
saveRDS(KEEP_FEATS, "data/interim/keep_features.rds")
cat("\n", paste(rule, collapse = "\n"), "\n")

#21

chk_cols <- c("AVG", "Tot size", "IPv", "LLC", "ARP")
t0 <- Sys.time()
d1 <- d2 <- d3 <- 0

for (i in seq_along(all_csv)) {
  d <- fread(all_csv[i], select = chk_cols, showProgress = FALSE)
  d1 <- max(d1, max(abs(d[["AVG"]] - d[["Tot size"]])))
  d2 <- max(d2, max(abs(d[["IPv"]] - d[["LLC"]])))
  d3 <- max(d3, max(abs(d[["ARP"]] - (1 - d[["IPv"]]))))
  rm(d)
  if (i %% 50 == 0 || i == length(all_csv))
    cat(sprintf("[%d/%d] %s\n", i, length(all_csv),
                format(round(difftime(Sys.time(), t0)))))
}

cat("\nacross all 46,776,697 rows:\n")
cat("max |AVG - Tot size| :", d1, "\n")
cat("max |IPv - LLC|      :", d2, "\n")
cat("max |ARP - (1-IPv)|  :", d3, "\n")

#22
library(ggplot2); library(scales)

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
  path <- file.path("figures", paste0(name, ".png"))
  ggsave(path, p, width = w, height = h, dpi = 200, bg = "white")
  cat("  ->", path, "\n"); invisible(path)
}

N_FULL      <- sum(lab_dt$rows)
BENIGN_FULL <- lab_dt[category == "Benign", sum(rows)]

# ---- Fig 1: all 34 labels ----
p1 <- ggplot(lab_dt, aes(reorder(label, rows), rows, fill = category)) +
  geom_col() + coord_flip() +
  scale_y_log10(labels = comma, expand = expansion(mult = c(0, 0.05))) +
  scale_fill_manual(values = PAL_8) +
  labs(title = "Fig 1. Distribution of all 34 CIC-IoT2023 labels",
       subtitle = "Log scale - four orders of magnitude separate the largest and smallest classes",
       x = NULL, y = "rows (log10)", fill = "Category",
       caption = sprintf("n = %s rows across %d files in %d label folders",
                         comma(N_FULL), sum(lab_dt$n_files), nrow(lab_dt))) +
  theme_ids()
save_fig(p1, "01_label_distribution", h = 8)

# ---- Fig 2: 8 classes ----
cat_plot <- copy(cat_dt)[, share := rows / sum(rows)]
p2 <- ggplot(cat_plot, aes(reorder(category, rows), rows, fill = category)) +
  geom_col(show.legend = FALSE) +
  geom_text(aes(label = sprintf("%s  (%.3f%%)", comma(rows), 100 * share)),
            hjust = -0.05, size = 3.2) +
  coord_flip() +
  scale_y_log10(labels = comma, expand = expansion(mult = c(0, 0.4))) +
  scale_fill_manual(values = PAL_8) +
  labs(title = "Fig 2. Rows per class (8-class taxonomy)",
       subtitle = "Web and BruteForce together are under 0.1% of the dataset",
       x = NULL, y = "rows (log10)") +
  theme_ids()
save_fig(p2, "02_category_distribution")

# ---- Fig 3: binary balance ----
bin_dt <- data.table(class = c("Benign", "Attack"),
                     count = c(BENIGN_FULL, N_FULL - BENIGN_FULL))
bin_dt[, share := count / sum(count)]
p3 <- ggplot(bin_dt, aes(x = "", y = share, fill = class)) +
  geom_col(width = 1) + coord_polar("y") +
  geom_text(aes(label = sprintf("%s\n%s\n%.4f%%", class, comma(count), 100 * share)),
            position = position_stack(vjust = 0.5), size = 4) +
  scale_fill_manual(values = c(Benign = "#2E7D32", Attack = "#C62828")) +
  labs(title = "Fig 3. Binary class balance",
       subtitle = "Why aggregate accuracy is uninformative on this dataset") +
  theme_void(base_size = 11) +
  theme(legend.position = "none",
        plot.title = element_text(face = "bold"),
        plot.subtitle = element_text(colour = "grey35"))
save_fig(p3, "03_binary_balance", w = 6, h = 6)

# ---- Fig 4: rows per source file ----
rc <- merge(rowcounts, lab_dt[, .(label, category)], by = "label")
p4 <- ggplot(rc, aes(rows, fill = category)) +
  geom_histogram(bins = 40, colour = "white", linewidth = 0.2) +
  scale_x_log10(labels = comma) +
  scale_fill_manual(values = PAL_8) +
  labs(title = "Fig 4. Rows per source file",
       subtitle = sprintf("%d files; sizes vary widely, so sampling quotas are allocated proportional to each file's row count",
                          nrow(rc)),
       x = "rows in file (log10)", y = "number of files", fill = "Category") +
  theme_ids()
save_fig(p4, "04_rows_per_file")

cat("\nfigures written to figures/\n")
print(list.files("figures"))

#23

# ---- Fig 5: sampling fidelity, attack half only ----
full_att <- cat_dt[category != "Benign", .(category, rows)]
full_att[, p_full := rows / sum(rows)]
ws_att <- ws[label_8 != "Benign", .N, by = label_8]
ws_att[, p_ws := N / sum(N)]
fid <- merge(full_att, ws_att, by.x = "category", by.y = "label_8")
fid[, pct_kept := 100 * N / rows]

fl <- melt(fid[, .(category, `Full dataset` = p_full, `Working set` = p_ws)],
           id.vars = "category", variable.name = "source", value.name = "share")

p5 <- ggplot(fl, aes(reorder(category, share), share, fill = source)) +
  geom_col(position = "dodge") + coord_flip() +
  scale_y_continuous(labels = percent) +
  scale_fill_manual(values = c(`Full dataset` = "grey55", `Working set` = "#1565C0")) +
  labs(title = "Fig 5. Attack-category shares: full dataset vs working set",
       subtitle = "A deliberate departure from proportional sampling - the rare categories are retained in full\nso that the composition arms have usable pools. The working set is not a miniature of the dataset.",
       x = NULL, y = "share of attack rows", fill = NULL) +
  theme_ids()
save_fig(p5, "05_sampling_fidelity")
print(fid[order(-rows), .(category, rows, kept = N,
                          pct_kept = round(pct_kept, 2),
                          p_full = round(100 * p_full, 2),
                          p_ws = round(100 * p_ws, 2))])

# ---- Fig 6: working-set composition ----
ws_comp <- ws[, .N, by = label_8]
p6 <- ggplot(ws_comp, aes(reorder(label_8, N), N, fill = label_8)) +
  geom_col(show.legend = FALSE) +
  geom_text(aes(label = comma(N)), hjust = -0.08, size = 3.2) +
  coord_flip() +
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.18))) +
  scale_fill_manual(values = PAL_8) +
  labs(title = "Fig 6. Working-set composition",
       subtitle = sprintf("%s rows x %d retained features, %.2f GB in memory; benign is %.1f%% by design",
                          comma(nrow(ws)), length(KEEP_FEATS),
                          as.numeric(object.size(ws)) / 1e9,
                          100 * mean(ws$label_binary == 0)),
       x = NULL, y = "rows") +
  theme_ids()
save_fig(p6, "06_working_set_composition")

# ---- Fig 7: feature constancy and what was dropped ----
cs <- copy(col_stats)
cs[, status := ifelse(column %in% DROP_COLS, "dropped (redundant)", "retained")]
p7 <- ggplot(cs, aes(reorder(column, modal_share), modal_share, fill = status)) +
  geom_col() + coord_flip() +
  geom_hline(yintercept = 0.999, linetype = "dashed", colour = "#C62828") +
  annotate("text", x = 3, y = 0.985, label = "0.999 threshold  ",
           hjust = 1, size = 3, colour = "#C62828") +
  scale_y_continuous(labels = percent, limits = c(0, 1.05)) +
  scale_fill_manual(values = c(`retained` = "grey60",
                               `dropped (redundant)` = "#C62828")) +
  labs(title = "Fig 7. Modal-value share per feature",
       subtitle = "No column has zero variance. The three dropped columns were removed as exact duplicates,\nnot for near-constancy - see Fig 8 for why the constancy threshold was not applied.",
       x = NULL, y = "share of rows holding the modal value", fill = NULL) +
  theme_ids()
save_fig(p7, "07_feature_constancy", h = 9)

# ---- Fig 8: why near-constant columns were kept ----
p8 <- ggplot(res, aes(label_8, reorder(column, nz_rate), fill = 100 * nz_rate)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%.2f", 100 * nz_rate)), size = 3, colour = "grey20") +
  scale_fill_gradient(low = "white", high = "#C62828", name = "% non-zero") +
  labs(title = "Fig 8. Sparse features are class-conditional, not uninformative",
       subtitle = "SSH fires in 32.1% of BruteForce rows against <=2.6% elsewhere. Removing near-constant\ncolumns would delete signal whose visibility depends on attack composition.",
       x = NULL, y = NULL) +
  theme_ids() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        panel.grid = element_blank())
save_fig(p8, "08_class_conditional_sparsity", w = 10, h = 5.5)

cat("\n"); print(list.files("figures"))

#24

fwrite(fid[order(-rows)], "reports/04_sampling_ledger.csv")
print(fid[order(-rows), .(category, full_rows = rows, kept = N,
                          pct_kept = round(100 * N / rows, 3),
                          p_full_pct = round(100 * p_full, 2),
                          p_ws_pct   = round(100 * p_ws, 2))])

set.seed(SEED)
sub <- ws[sample.int(.N, 200000)]

# ---- Fig 9: correlation among retained features ----
cmk  <- cm[KEEP_FEATS, KEEP_FEATS]
ord  <- hclust(as.dist(1 - abs(cmk)))$order
lv   <- KEEP_FEATS[ord]
cdt2 <- as.data.table(as.table(cmk)); setnames(cdt2, c("f1", "f2", "r"))
cdt2[, f1 := factor(as.character(f1), levels = lv)]
cdt2[, f2 := factor(as.character(f2), levels = lv)]

p9 <- ggplot(cdt2, aes(f1, f2, fill = r)) + geom_tile() +
  scale_fill_gradient2(low = "#1565C0", mid = "white", high = "#C62828",
                       midpoint = 0, limits = c(-1, 1)) +
  labs(title = "Fig 9. Correlation among the 36 retained features",
       subtitle = "Clustered ordering. Redundant blocks motivate selection; exact duplicates were already removed.",
       x = NULL, y = NULL, fill = "Pearson r") +
  theme_ids() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
        axis.text.y = element_text(size = 7),
        panel.grid = element_blank())
save_fig(p9, "09_correlation_heatmap", w = 10, h = 9)

skew1 <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 3) return(0)
  s <- sd(x); if (s == 0) return(0)
  mean((x - mean(x))^3) / s^3
}
sk <- data.table(column = KEEP_FEATS,
                 skew = sapply(KEEP_FEATS, function(c) skew1(ws[[c]])))
sk <- sk[order(-abs(skew))]          # <- was setorder(sk, -abs(skew))
print(sk[1:12])
fwrite(sk, "reports/04_skewness.csv")

p10 <- ggplot(sk, aes(reorder(column, abs(skew)), skew, fill = abs(skew) > 5)) +
  geom_col() + coord_flip() +
  scale_fill_manual(values = c(`FALSE` = "grey60", `TRUE` = "#C62828"),
                    labels = c("moderate", "heavy-tailed"), name = NULL) +
  labs(title = "Fig 10. Skewness per retained feature",
       subtitle = "Heavy tails are why scaling uses median/MAD rather than mean/SD, with clipping at +/-3",
       x = NULL, y = "skewness") +
  theme_ids()
save_fig(p10, "10_skewness", h = 9)

set.seed(SEED)
sub <- ws[sample.int(.N, 200000)]
top <- sk[1:12, column]
long <- melt(sub[, c(top, "label_binary"), with = FALSE],
             id.vars = "label_binary", variable.name = "feature",
             value.name = "value")
long <- long[is.finite(value)]
long[, cls := ifelse(label_binary == 1, "Attack", "Benign")]

p11 <- ggplot(long, aes(value, fill = cls)) +
  geom_density(alpha = 0.5, colour = NA) +
  facet_wrap(~feature, scales = "free", ncol = 4) +
  scale_fill_manual(values = c(Benign = "#2E7D32", Attack = "#C62828")) +
  labs(title = "Fig 11. Feature distributions before scaling",
       subtitle = "The twelve most heavy-tailed retained features, benign vs attack",
       x = NULL, y = "density", fill = NULL) +
  theme_ids(9)
save_fig(p11, "11_distributions_raw", w = 12, h = 8)

# the sampling ledger, printed properly this time
print(fid[order(-rows), .(category, full_rows = rows, kept = N,
                          pct_kept = round(100 * N / rows, 3),
                          p_full_pct = round(100 * p_full, 2),
                          p_ws_pct   = round(100 * p_ws, 2))])

