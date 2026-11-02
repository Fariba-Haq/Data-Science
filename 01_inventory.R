


RAW_DIR <- "D:/0.1.Researches/DS/1.Dataset_author"
SEED    <- 13

OUT_DIR <- "data"
FIG_DIR <- "figures"
REP_DIR <- "reports"


for (p in c("data.table", "ggplot2", "scales"))
  if (!requireNamespace(p, quietly = TRUE))
    stop(sprintf("Install it first:  install.packages('%s')", p))

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(scales)
})

for (d in c(OUT_DIR, file.path(OUT_DIR, "interim"), file.path(OUT_DIR, "processed"),
            FIG_DIR, REP_DIR))
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

CAT_MAP <- c(
  "DDoS-ICMP_Flood"="DDoS", "DDoS-UDP_Flood"="DDoS", "DDoS-TCP_Flood"="DDoS",
  "DDoS-PSHACK_FLOOD"="DDoS", "DDoS-SYN_Flood"="DDoS", "DDoS-RSTFINFLOOD"="DDoS",
  "DDoS-SynonymousIP_Flood"="DDoS", "DDoS-ICMP_Fragmentation"="DDoS",
  "DDoS-UDP_Fragmentation"="DDoS", "DDoS-ACK_Fragmentation"="DDoS",
  "DDoS-HTTP_Flood"="DDoS", "DDoS-SlowLoris"="DDoS",
  "DoS-UDP_Flood"="DoS", "DoS-TCP_Flood"="DoS",
  "DoS-SYN_Flood"="DoS", "DoS-HTTP_Flood"="DoS",
  "Mirai-greeth_flood"="Mirai", "Mirai-greip_flood"="Mirai",
  "Mirai-udpplain"="Mirai",
  "VulnerabilityScan"="Recon", "Recon-HostDiscovery"="Recon",
  "Recon-OSScan"="Recon", "Recon-PortScan"="Recon", "Recon-PingSweep"="Recon",
  "MITM-ArpSpoofing"="Spoofing", "DNS_Spoofing"="Spoofing",
  "BrowserHijacking"="Web", "CommandInjection"="Web", "SqlInjection"="Web",
  "XSS"="Web", "Backdoor_Malware"="Web", "Uploading_Attack"="Web",
  "DictionaryBruteForce"="BruteForce",
  "Benign_Final"="Benign"
)
BENIGN_FOLDER <- "Benign_Final"
ATTACK_CATS   <- c("DDoS","DoS","Mirai","Recon","Spoofing","Web","BruteForce")


if (!dir.exists(RAW_DIR)) stop("RAW_DIR not found: ", RAW_DIR)

all_csv <- list.files(RAW_DIR, pattern = "\\.csv$", full.names = TRUE,
                      recursive = TRUE, ignore.case = TRUE)
if (!length(all_csv)) stop("No .csv files found under ", RAW_DIR)

cat(sprintf("Stage 1: %d CSV files in %d folders\n",
            length(all_csv), uniqueN(basename(dirname(all_csv)))))


hdrs <- lapply(all_csv, function(f) names(fread(f, nrows = 0)))
sig  <- vapply(hdrs, paste, "", collapse = "|")
u    <- unique(sig)

cat("distinct header signatures:", length(u), "\n")
if (length(u) != 1L) {
  print(table(match(sig, u)))
  stop("Files do not share one schema. Resolve before continuing.")
}

ref <- hdrs[[1]]  
cat("columns:", length(ref), "\n")
cat("label column present:",
    any(tolower(ref) %in% c("label", "class")), "\n")


t0  <- Sys.time()
res <- vector("list", length(all_csv))

for (i in seq_along(all_csv)) {
  n <- tryCatch(nrow(fread(all_csv[i], select = 1L, showProgress = FALSE)),
                error = function(e) NA_integer_)
  res[[i]] <- data.table(label = basename(dirname(all_csv[i])),
                         file  = basename(all_csv[i]),
                         rows  = n)
  if (i %% 50 == 0 || i == length(all_csv))
    cat(sprintf("  [%d/%d] %s\n", i, length(all_csv),
                format(round(difftime(Sys.time(), t0)))))
}

rowcounts <- rbindlist(res)
if (any(is.na(rowcounts$rows))) stop("Some files could not be read.")



lab_dt <- rowcounts[, .(n_files = .N, rows = sum(rows)), by = label]
lab_dt[, category  := CAT_MAP[label]]
lab_dt[, share_pct := round(100 * rows / sum(rows), 4)]
setorder(lab_dt, -rows)

miss_map <- setdiff(lab_dt$label, names(CAT_MAP))
miss_dir <- setdiff(names(CAT_MAP), lab_dt$label)
if (length(miss_map)) stop("Folders with no CAT_MAP entry: ", paste(miss_map, collapse = ", "))
if (length(miss_dir)) stop("CAT_MAP entries with no folder: ", paste(miss_dir, collapse = ", "))
cat("taxonomy verified: all", nrow(lab_dt), "folders mapped\n")

cat_dt <- lab_dt[, .(n_labels = .N, n_files = sum(n_files), rows = sum(rows)),
                 by = category]
cat_dt[, share_pct := round(100 * rows / sum(rows), 4)]
setorder(cat_dt, -rows)

N_FULL      <- sum(lab_dt$rows)
BENIGN_FULL <- lab_dt[category == "Benign", sum(rows)]


summary_txt <- c(
  "CIC-IoT2023 — STAGE 1 INVENTORY",
  sprintf("run at        : %s", format(Sys.time())),
  sprintf("R / data.table: %s / %s", getRversion(),
          as.character(packageVersion("data.table"))),
  sprintf("RAW_DIR       : %s", normalizePath(RAW_DIR)),
  "",
  sprintf("csv files     : %d", nrow(rowcounts)),
  sprintf("label folders : %d", nrow(lab_dt)),
  sprintf("total rows    : %s", format(N_FULL, big.mark = ",")),
  sprintf("raw columns   : %d (no label column; label = folder name)", length(ref)),
  sprintf("benign rows   : %s", format(BENIGN_FULL, big.mark = ",")),
  sprintf("BENIGN SHARE  : %.4f%%", 100 * BENIGN_FULL / N_FULL),
  sprintf("smallest label: %s = %s rows",
          lab_dt[which.min(rows), label], format(min(lab_dt$rows), big.mark = ",")),
  "",
  "--- columns ---",
  paste(sprintf("%2d. %s", seq_along(ref), ref), collapse = "\n"),
  "",
  "--- categories ---",
  paste(capture.output(print(cat_dt)), collapse = "\n")
)
writeLines(summary_txt, file.path(REP_DIR, "01_inventory_summary.txt"))
cat("\n", paste(summary_txt[1:14], collapse = "\n"), "\n\n", sep = "")

fwrite(rowcounts, file.path(REP_DIR, "01_rows_per_file.csv"))
fwrite(lab_dt,    file.path(REP_DIR, "01_label_counts.csv"))
fwrite(cat_dt,    file.path(REP_DIR, "01_category_counts.csv"))


## Fig 1

p1 <- ggplot(lab_dt, aes(reorder(label, rows), rows, fill = category)) +
  geom_col() + coord_flip() +
  scale_y_log10(labels = comma, expand = expansion(mult = c(0, 0.05))) +
  scale_fill_manual(values = PAL_8) +
  labs(title = "Fig 1. Distribution of all 34 CIC-IoT2023 labels",
       subtitle = "Log scale - four orders of magnitude separate the largest and smallest classes",
       x = NULL, y = "rows (log10)", fill = "Category",
       caption = sprintf("n = %s rows across %d files in %d label folders",
                         comma(N_FULL), nrow(rowcounts), nrow(lab_dt))) +
  theme_ids()
save_fig(p1, "01_label_distribution", h = 8)

## Fig 2
cat_plot <- copy(cat_dt)[, share := rows / sum(rows)]
p2 <- ggplot(cat_plot, aes(reorder(category, rows), rows, fill = category)) +
  geom_col(show.legend = FALSE) +
  geom_text(aes(label = sprintf("%s  (%.3f%%)", comma(rows), 100 * share)),
            hjust = -0.05, size = 3.2) +
  coord_flip() +
  scale_y_log10(labels = comma, expand = expansion(mult = c(0, 0.15))) +
  scale_fill_manual(values = PAL_8) +
  labs(title = "Fig 2. Rows per class (8-class taxonomy)",
       subtitle = "Web and BruteForce together are under 0.1% of the dataset",
       x = NULL, y = "rows (log10)") +
  theme_ids()
save_fig(p2, "02_category_distribution")

## Fig 3
bin_dt <- data.table(class = c("Benign", "Attack"),
                     count = c(BENIGN_FULL, N_FULL - BENIGN_FULL))
bin_dt[, share := count / sum(count)]
p3 <- ggplot(bin_dt, aes(x = "", y = share, fill = class)) +
  geom_col(width = 1) + coord_polar("y") +
  geom_text(aes(label = sprintf("%s\n%s\n%.4f%%", class, comma(count), 100 * share)),
            position = position_stack(vjust = 0.5), size = 4) +
  scale_fill_manual(values = c(Benign = "#2E7D32", Attack = "#C62828")) +
  labs(title = "Fig 3. Binary class balance",
       subtitle = "Predicting 'attack' for every row scores 97.65% accuracy - which is why we report macro-F1") +
  theme_void(base_size = 11) +
  theme(legend.position = "none",
        plot.title = element_text(face = "bold"),
        plot.subtitle = element_text(colour = "grey35"))
save_fig(p3, "03_binary_balance", w = 6, h = 6)

## Fig 4
rc <- merge(rowcounts, lab_dt[, .(label, category)], by = "label")
p4 <- ggplot(rc, aes(rows, fill = category)) +
  geom_histogram(bins = 40, colour = "white", linewidth = 0.2) +
  scale_x_log10(labels = comma) +
  scale_fill_manual(values = PAL_8) +
  labs(title = "Fig 4. Rows per source file",
       subtitle = sprintf("%d files spanning three orders of magnitude, so sampling quotas are allocated proportional to each file's row count",
                          nrow(rc)),
       x = "rows in file (log10)", y = "number of files", fill = "Category") +
  theme_ids()
save_fig(p4, "04_rows_per_file")



saveRDS(list(all_csv = all_csv, ref = ref, rowcounts = rowcounts,
             lab_dt = lab_dt, cat_dt = cat_dt,
             N_FULL = N_FULL, BENIGN_FULL = BENIGN_FULL,
             CAT_MAP = CAT_MAP, RAW_DIR = RAW_DIR),
        file.path(OUT_DIR, "interim", "inventory.rds"))

cat("\nStage 1 complete. -> data/interim/inventory.rds\n")
