

SEED            <- 13
CONFIG_ROWS     <- 30000     
BENIGN_FRACTION <- 0.5        
BASE_MIX        <- "uniform" 

N_CONTROL       <- 20      
N_SEEDS_DROP    <- 10        
N_MIXES         <- 25       
N_SEEDS_DIR     <- 5         
DIRICHLET_CONC  <- 0.4    

OUT_DIR <- "data"; FIG_DIR <- "figures"; REP_DIR <- "reports"

ATTACK_CATS <- c("DDoS","DoS","Mirai","Recon","Spoofing","Web","BruteForce")



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
PAL_8 <- c(Benign = "#2E7D32", DDoS = "#C62828", DoS = "#EF6C00",
           Mirai = "#6A1B9A", Recon = "#1565C0", Spoofing = "#00838F",
           Web = "#AD1457", BruteForce = "#5D4037")
save_fig <- function(p, name, w = 9, h = 5.5) {
  path <- file.path(FIG_DIR, paste0(name, ".png"))
  ggsave(path, p, width = w, height = h, dpi = 200, bg = "white")
  cat("  figure ->", path, "\n"); invisible(path)
}


sp  <- readRDS(file.path(OUT_DIR, "processed", "split_safe.rds"))
inv <- readRDS(file.path(OUT_DIR, "interim", "inventory.rds"))
y8  <- sp$y8_train

pool <- lapply(c("Benign", ATTACK_CATS), function(cl) which(y8 == cl))
names(pool) <- c("Benign", ATTACK_CATS)
pool_n <- sapply(pool, length)

n_benign <- round(CONFIG_ROWS * BENIGN_FRACTION)
n_attack <- CONFIG_ROWS - n_benign

cat(sprintf("Stage 6: %s training rows | config = %s (%s benign + %s attack)\n",
            format(length(y8), big.mark = ","), comma(CONFIG_ROWS),
            comma(n_benign), comma(n_attack)))
print(data.table(class = names(pool_n), pool = as.integer(pool_n),
                 max_share_of_attack_half =
                   round(pmin(1, pool_n / n_attack), 3)))

if (pool_n[["Benign"]] < n_benign)
  stop("Benign pool too small. Lower CONFIG_ROWS or BENIGN_FRACTION.")


cat_full <- inv$cat_dt[category %in% ATTACK_CATS]
natural  <- setNames(numeric(length(ATTACK_CATS)), ATTACK_CATS)
natural[cat_full$category] <- cat_full$rows
natural <- natural / sum(natural)

uniform <- setNames(rep(1 / length(ATTACK_CATS), length(ATTACK_CATS)), ATTACK_CATS)

mix_dt <- data.table(category = ATTACK_CATS,
                     natural_pct = round(100 * natural[ATTACK_CATS], 4),
                     uniform_pct = round(100 * uniform[ATTACK_CATS], 4),
                     rows_if_natural = round(natural[ATTACK_CATS] * n_attack),
                     rows_if_uniform = round(uniform[ATTACK_CATS] * n_attack))
cat("\n--- candidate base mixes, rows per category in the attack half ---\n")
print(mix_dt)
fwrite(mix_dt, file.path(REP_DIR, "06_base_mixes.csv"))

base_w <- switch(BASE_MIX, uniform = uniform, natural = natural,
                 stop("BASE_MIX must be 'uniform' or 'natural'"))
cat(sprintf("\nBASE_MIX = %s\n", BASE_MIX))


js_divergence <- function(p, q) {
  p <- p / sum(p); q <- q / sum(q); m <- 0.5 * (p + q)
  kl <- function(a, b) { i <- a > 0; sum(a[i] * log2(a[i] / b[i])) }
  0.5 * kl(p, m) + 0.5 * kl(q, m)
}

allocate <- function(w, n_attack, pool_n) {
  w <- w / sum(w)
  q <- floor(w * n_attack)
  short <- n_attack - sum(q)                    
  if (short > 0) {
    ord <- order(w * n_attack - q, decreasing = TRUE)
    q[ord[seq_len(short)]] <- q[ord[seq_len(short)]] + 1
  }
  cap <- pool_n[names(q)]
  repeat {
    over <- q > cap
    if (!any(over)) break
    deficit <- sum(q[over] - cap[over])
    q[over]  <- cap[over]
    room <- cap - q
    can  <- room > 0 & !over
    if (!any(can) || deficit == 0) break
    wts <- w[can] / sum(w[can])
    add <- pmin(room[can], floor(wts * deficit))
    q[can] <- q[can] + add
    left <- deficit - sum(add)
    if (left > 0) {                                
      for (nm in names(sort(room[can] - add, decreasing = TRUE))) {
        if (left == 0) break
        give <- min(left, room[[nm]] - add[[nm]])
        q[[nm]] <- q[[nm]] + give; left <- left - give
      }
    }
    if (sum(q) >= n_attack) break
  }
  q
}

draw_config <- function(w, seed) {
  set.seed(seed)
  q   <- allocate(w, n_attack, pool_n[ATTACK_CATS])
  idx <- sample(pool[["Benign"]], n_benign)
  for (cl in ATTACK_CATS)
    if (q[[cl]] > 0) idx <- c(idx, sample(pool[[cl]], q[[cl]]))
  list(idx = sample(idx), q = q)                   # shuffled
}




specs <- list()

for (b in seq_len(N_CONTROL))
  specs[[length(specs) + 1]] <- list(arm = "control",
                                     config = "control_base",
                                     w = base_w, seed = SEED + b)

for (cl in ATTACK_CATS) {
  w <- base_w; w[cl] <- 0; w <- w / sum(w)
  for (s in seq_len(N_SEEDS_DROP))
    specs[[length(specs) + 1]] <- list(arm = "drop_one",
                                       config = paste0("drop_", cl),
                                       w = w, seed = SEED + 100 * s + which(ATTACK_CATS == cl))
}

set.seed(SEED)
alpha <- pmax(base_w * length(ATTACK_CATS) * DIRICHLET_CONC, 0.05)
dir_w <- lapply(seq_len(N_MIXES), function(m) {
  g <- rgamma(length(alpha), shape = alpha, rate = 1)
  setNames(g / sum(g), ATTACK_CATS)
})
for (m in seq_len(N_MIXES))
  for (s in seq_len(N_SEEDS_DIR))
    specs[[length(specs) + 1]] <- list(arm = "dirichlet",
                                       config = sprintf("dirichlet_%02d", m),
                                       w = dir_w[[m]], seed = SEED + 1000 * s + m)


for (s in seq_len(N_SEEDS_DROP))
  specs[[length(specs) + 1]] <- list(arm = "reference",
                                     config = "natural_mix",
                                     w = natural, seed = SEED + 5000 + s)

cat(sprintf("\nbuilding %d configurations ...\n", length(specs)))

draws <- vector("list", length(specs))
diag  <- vector("list", length(specs))
t0 <- Sys.time()

for (i in seq_along(specs)) {
  s  <- specs[[i]]
  dc <- draw_config(s$w, s$seed)
  key <- sprintf("%s__seed%d", s$config, s$seed)
  draws[[i]] <- dc$idx
  names(draws)[i] <- key

  qv       <- as.numeric(dc$q[ATTACK_CATS])
  realised <- setNames(qv / max(sum(qv), 1), ATTACK_CATS)
  target   <- setNames(as.numeric(s$w[ATTACK_CATS]), ATTACK_CATS)
  short    <- ATTACK_CATS[qv < floor(target * n_attack) - 1]

  diag[[i]] <- cbind(
    data.table(
      draw_key = key, arm = s$arm, config = s$config, seed = s$seed,
      n = length(dc$idx),
      benign_frac = n_benign / length(dc$idx),
      js_target_vs_base   = js_divergence(target,   base_w[ATTACK_CATS]),
      js_realised_vs_base = js_divergence(realised, base_w[ATTACK_CATS]),
      shortfall = paste(short, collapse = ";")
    ),
    as.data.table(as.list(realised))
  )
  if (i %% 50 == 0 || i == length(specs))
    cat(sprintf("  [%d/%d] %s\n", i, length(specs),
                format(round(difftime(Sys.time(), t0)))))
}

diag <- rbindlist(diag)


ok_n  <- uniqueN(diag$n) == 1
ok_bf <- uniqueN(round(diag$benign_frac, 9)) == 1

cat("\n--- INVARIANTS ---\n")
cat("n identical across all configurations      :", ok_n,
    sprintf("(%s)\n", paste(unique(diag$n), collapse = ", ")))
cat("benign fraction identical across all configs:", ok_bf,
    sprintf("(%.6f)\n", diag$benign_frac[1]))
if (!ok_n || !ok_bf)
  stop("Invariant violated: composition would be confounded with n or imbalance.")

n_short <- sum(nzchar(diag$shortfall))
cat("configurations with a quota shortfall      :", n_short,
    sprintf("(%.1f%%)\n", 100 * n_short / nrow(diag)))
if (n_short) {
  cat("  categories involved:",
      paste(sort(unique(unlist(strsplit(diag$shortfall[nzchar(diag$shortfall)], ";")))),
            collapse = ", "), "\n")
  cat("  -> report REALISED weights, not target weights, for these.\n")
}

cat("\n--- per-arm summary ---\n")
print(diag[, .(configs = .N,
               js_mean = round(mean(js_realised_vs_base), 4),
               js_max  = round(max(js_realised_vs_base), 4),
               shortfalls = sum(nzchar(shortfall))), by = arm])

fwrite(diag, file.path(REP_DIR, "06_arm_diagnostics.csv"))

saveRDS(list(draws = draws, diag = diag,
             base_mix = BASE_MIX, base_w = base_w,
             natural_w = natural, uniform_w = uniform,
             config_rows = CONFIG_ROWS, benign_fraction = BENIGN_FRACTION,
             n_benign = n_benign, n_attack = n_attack, seed = SEED,
             note = paste("draws are row indices into X_train / y8_train of",
                          "data/processed/split_safe.rds")),
        file.path(OUT_DIR, "processed", "arms.rds"))
cat("\nsaved -> data/processed/arms.rds\n")


##  Fig 19: realised composition of each configuration
first_seed <- diag[, .SD[1], by = config]
cl <- melt(first_seed[, c("config", "arm", ATTACK_CATS), with = FALSE],
           id.vars = c("config", "arm"), variable.name = "category",
           value.name = "share")
cl[, arm := factor(arm, levels = c("control", "reference", "drop_one", "dirichlet"))]

p19 <- ggplot(cl[arm %in% c("control", "reference", "drop_one")],
              aes(config, share, fill = category)) +
  geom_col() + coord_flip() +
  scale_y_continuous(labels = percent) +
  scale_fill_manual(values = PAL_8[ATTACK_CATS]) +
  labs(title = "Fig 19. Attack-mix composition per configuration",
       subtitle = sprintf("Total rows (%s) and benign:attack ratio (%.0f:%.0f) are IDENTICAL across every configuration.\nOnly the internal composition of the attack half varies.",
                          comma(CONFIG_ROWS), 100 * BENIGN_FRACTION,
                          100 * (1 - BENIGN_FRACTION)),
       x = NULL, y = "share of the attack half", fill = NULL) +
  theme_ids()
save_fig(p19, "19_composition_arms", h = 6)

##  Fig 20: the dose-response x-axis
p20 <- ggplot(diag[arm == "dirichlet"], aes(js_realised_vs_base)) +
  geom_histogram(bins = 20, fill = "#6A1B9A", colour = "white") +
  labs(title = "Fig 20. Divergence of the sampled mixes from the base composition",
       subtitle = "This is the x-axis of the dose-response curve in Results. A spread of values is what\nmakes the relationship estimable; a cluster would leave the curve unidentified.",
       x = "Jensen-Shannon divergence from the base mix",
       y = "configurations") +
  theme_ids()
save_fig(p20, "20_js_divergence")

## Fig 21: the design, in one picture
p21 <- ggplot(diag, aes(reorder(arm, js_realised_vs_base), js_realised_vs_base,
                        colour = arm)) +
  geom_jitter(width = 0.18, alpha = 0.6, size = 1.4, show.legend = FALSE) +
  stat_summary(fun = median, geom = "crossbar", width = 0.5,
               colour = "grey20", linewidth = 0.3, show.legend = FALSE) +
  coord_flip() +
  labs(title = "Fig 21. Perturbation magnitude by arm",
       subtitle = "The control arm sits at zero divergence BY CONSTRUCTION - it varies the random draw, not the mix.\nAny subset movement it shows is sampling noise, and is the yardstick the other arms are read against.",
       x = NULL, y = "Jensen-Shannon divergence from the base mix") +
  theme_ids()
save_fig(p21, "21_perturbation_magnitude", h = 4.5)

si <- c(
  sprintf("timestamp       : %s", format(Sys.time())),
  sprintf("R version       : %s", R.version.string),
  sprintf("seed            : %d", SEED),
  sprintf("BASE_MIX        : %s", BASE_MIX),
  sprintf("CONFIG_ROWS     : %d  (%d benign + %d attack)", CONFIG_ROWS, n_benign, n_attack),
  sprintf("configurations  : %d", nrow(diag)),
  sprintf("  control       : %d", sum(diag$arm == "control")),
  sprintf("  drop_one      : %d", sum(diag$arm == "drop_one")),
  sprintf("  dirichlet     : %d", sum(diag$arm == "dirichlet")),
  sprintf("  reference     : %d", sum(diag$arm == "reference")),
  sprintf("n identical     : %s", ok_n),
  sprintf("benign frac id. : %s", ok_bf),
  sprintf("shortfalls      : %d", n_short),
  sprintf("max JS from base: %.4f", max(diag$js_realised_vs_base))
)
writeLines(si, file.path(REP_DIR, "06_run_record.txt"))
cat("\n", paste(si, collapse = "\n"), "\n")
cat("\nStage 6 complete. -> data/processed/arms.rds\n")
