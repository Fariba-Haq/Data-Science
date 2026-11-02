


RAW_DIR <- "D:/0.1.Researches/DS/1.Dataset_author"
SEED    <- 13

OUT_DIR <- "data"
FIG_DIR <- "figures"
REP_DIR <- "reports"

CATEGORY_CAP    <- 200000   
MIN_PER_LABEL   <- 3000     
BENIGN_TAKE_ALL <- TRUE     



suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(scales)
})

for (d in c(OUT_DIR, file.path(OUT_DIR, "interim"), FIG_DIR, REP_DIR))
  dir.create(d, showWarnings = FALSE, recursive = TRUE)

theme_ids <- function(base = 11) {
  theme_minimal(base_size = base) +
    theme(panel.grid.minor = element_blank(),
          plot.title    = element_text(face = "bold", size = base + 2),
          plot.subtitle = element_text(colour = "grey35", size = base - 1),
          plot.caption  = element_text(colour = "grey45", size = base - 2, hjust = 0),
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

inv <- readRDS(file.path(OUT_DIR, "interim", "inventory.rds"))
all_csv   <- inv$all_csv
ref       <- inv$ref
rowcounts <- inv$rowcounts
lab_dt    <- copy(inv$lab_dt)
cat_dt    <- copy(inv$cat_dt)
CAT_MAP   <- inv$CAT_MAP
BENIGN_FOLDER <- "Benign_Final"

cat(sprintf("Stage 2: %s source rows across %d files\n",
            format(inv$N_FULL, big.mark = ","), length(all_csv)))



lab_dt[, cat_rows   := sum(rows), by = category]
lab_dt[, cat_target := pmin(cat_rows, CATEGORY_CAP)]
if (BENIGN_TAKE_ALL) lab_dt[category == "Benign", cat_target := cat_rows]

lab_dt[, quota := pmin(rows,
                       pmax(round(cat_target * rows / cat_rows),
                            pmin(MIN_PER_LABEL, rows)))]

WS_ROWS <- sum(lab_dt$quota)
cat(sprintf("planned working set: %s rows (~%.2f GB numeric)\n",
            format(WS_ROWS, big.mark = ","), WS_ROWS * length(ref) * 8 / 1e9))

print(lab_dt[order(category, -rows),
             .(label, category, rows, quota,
               pct_kept = round(100 * quota / rows, 2))], nrows = 40)


plan <- merge(rowcounts, lab_dt[, .(label, label_rows = rows, quota)], by = "label")
setorder(plan, label, -rows, file)            

plan[, take := round(quota * rows / label_rows)]
plan[, take := { t <- take; t[1] <- t[1] + (quota[1] - sum(t)); t }, by = label]
plan[, take := pmax(pmin(take, rows), 0L)]

chk <- plan[, .(allocated = sum(take), quota = quota[1]), by = label][allocated != quota]
if (nrow(chk)) { print(chk); stop("File allocation does not sum to quota.") }
cat("allocation verified; rows to read:", format(sum(plan$take), big.mark = ","), "\n")
fwrite(plan, file.path(REP_DIR, "02_sampling_plan.csv"))


set.seed(SEED)
t0    <- Sys.time()
parts <- vector("list", nrow(plan))

for (i in seq_len(nrow(plan))) {
  if (plan$take[i] <= 0) next
  f  <- file.path(RAW_DIR, plan$label[i], plan$file[i])
  dt <- fread(f, showProgress = FALSE)

  stopifnot(identical(names(dt), ref))          

  if (plan$take[i] < nrow(dt))
    dt <- dt[sample.int(nrow(dt), plan$take[i])]

  dt[, label_fine := plan$label[i]]             
  parts[[i]] <- dt
  rm(dt)

  if (i %% 50 == 0 || i == nrow(plan)) {
    gc(verbose = FALSE)
    cat(sprintf("  [%d/%d] %s\n", i, nrow(plan),
                format(round(difftime(Sys.time(), t0)))))
  }
}

ws <- rbindlist(parts, use.names = TRUE)
rm(parts); gc(verbose = FALSE)

ws[, label_8      := CAT_MAP[label_fine]]                   
ws[, label_binary := as.integer(label_fine != BENIGN_FOLDER)] 
ws <- ws[sample.int(.N)]     

saveRDS(ws, file.path(OUT_DIR, "interim", "working_set.rds"))


got <- ws[, .N, by = label_fine]
vfy <- merge(lab_dt[, .(label_fine = label, quota)], got, by = "label_fine")
vfy[, delta := N - quota]
if (any(vfy$delta != 0)) { print(vfy[delta != 0]); stop("Label counts off quota.") }
if (any(is.na(ws$label_8)))                            stop("Unmapped label_8 present.")
cat("verified: every label on quota, no unmapped categories\n")

cat(sprintf("rows %s | cols %d | %.2f GB | benign %.2f%%\n",
            format(nrow(ws), big.mark = ","), ncol(ws),
            as.numeric(object.size(ws)) / 1e9,
            100 * mean(ws$label_binary == 0)))


## Fig 5 
full_att <- cat_dt[category != "Benign", .(category, rows)]
full_att[, p_full := rows / sum(rows)]
ws_att   <- ws[label_8 != "Benign", .N, by = label_8]
ws_att[, p_ws := N / sum(N)]
fid <- merge(full_att, ws_att, by.x = "category", by.y = "label_8")

fl <- melt(fid[, .(category, `Full dataset` = p_full, `Working set` = p_ws)],
           id.vars = "category", variable.name = "source", value.name = "share")
p5 <- ggplot(fl, aes(reorder(category, share), share, fill = source)) +
  geom_col(position = "dodge") + coord_flip() +
  scale_y_continuous(labels = percent) +
  scale_fill_manual(values = c(`Full dataset` = "grey55", `Working set` = "#1565C0")) +
  labs(title = "Fig 5. Attack-category shares: full dataset vs working set",
       subtitle = "A deliberate departure from proportional sampling - rare categories are retained in full\nso the composition arms have usable pools. The working set is not a miniature of the dataset.",
       x = NULL, y = "share of attack rows", fill = NULL) +
  theme_ids()
save_fig(p5, "05_sampling_fidelity")

print(fid[order(-rows), .(category, full_rows = rows, kept = N,
                          pct_kept   = round(100 * N / rows, 3),
                          p_full_pct = round(100 * p_full, 2),
                          p_ws_pct   = round(100 * p_ws, 2))])
fwrite(fid, file.path(REP_DIR, "02_sampling_ledger.csv"))

## Fig 6 — what we actually carry forward.
p6 <- ggplot(ws[, .N, by = label_8], aes(reorder(label_8, N), N, fill = label_8)) +
  geom_col(show.legend = FALSE) +
  geom_text(aes(label = comma(N)), hjust = -0.08, size = 3.2) +
  coord_flip() +
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.18))) +
  scale_fill_manual(values = PAL_8) +
  labs(title = "Fig 6. Working-set composition",
       subtitle = sprintf("%s rows, %.2f GB in memory; benign is %.1f%% by design",
                          comma(nrow(ws)), as.numeric(object.size(ws)) / 1e9,
                          100 * mean(ws$label_binary == 0)),
       x = NULL, y = "rows") +
  theme_ids()
save_fig(p6, "06_working_set_composition")

si <- c(
  sprintf("timestamp      : %s", format(Sys.time())),
  sprintf("R version      : %s", R.version.string),
  sprintf("platform       : %s", R.version$platform),
  sprintf("data.table     : %s", as.character(packageVersion("data.table"))),
  sprintf("seed           : %d", SEED),
  sprintf("CATEGORY_CAP   : %d", CATEGORY_CAP),
  sprintf("MIN_PER_LABEL  : %d", MIN_PER_LABEL),
  sprintf("BENIGN_TAKE_ALL: %s", BENIGN_TAKE_ALL),
  sprintf("source rows    : %d", inv$N_FULL),
  sprintf("working set    : %d", nrow(ws))
)
writeLines(si, file.path(REP_DIR, "02_run_record.txt"))
cat("\n", paste(si, collapse = "\n"), "\n")
cat("\nStage 2 complete. -> data/interim/working_set.rds\n")
