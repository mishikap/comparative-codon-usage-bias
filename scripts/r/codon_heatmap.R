## ---- Required packages you already have ----
library(readr)    # read_lines()
library(stringr)  # str_detect()
library(lattice)  # levelplot()

## ========= USER SETTINGS =========

# 1) Folder that contains six TXT files
data_dir <- file.path("data", "cusp")  

# 2) Exact filenames (as they appear in that folder) + pretty labels
files <- c(
  "K12codonusage.txt",
  "Bsubtilis168codonusage.txt",
  "CodoncountSentericaLT2.txt",
  "CodoncountPputida.txt",
  "CodoncountHsalinariumNRC-1.txt",
  "CodoncountScerevisiaS288C.txt"
)

organism_labels <- c("E. coli K-12",
                     "B. subtilis 168",
                     "S. enterica LT2",
                     "P. putida KT2440",
                     "H. salinarium NRC-1",
                     "S. cerevisiae S288C")

# Toggle: use RSCU (TRUE) or per-AA normalized frequency (FALSE)
use_rscu <- TRUE
top_variance_N <- 40
every_n_codon_label <- 1

# Where to save figures (absolute path)
out_dir <- file.path("figures", "generated_r")

## ========= END USER SETTINGS =========

# Create output folder if needed
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# Helper: parse one codon-usage file in your format
# Accepts header "#Codon AA Fraction Frequency Number" after some comment/blank lines.
read_codon_table <- function(file_path, org_label){
  raw <- read_lines(file_path)

  # find header (allow with/without '#' in front)
  hdr_idx <- which(str_detect(raw, "^#?Codon\\s+AA\\s+Fraction\\s+Frequency\\s+Number"))[1]
  if (is.na(hdr_idx)) stop("Header not found in file: ", file_path)

  # table starts after header, ends at first blank or comment or EOF
  tab_lines <- raw[(hdr_idx+1):length(raw)]
  end_idx <- which(tab_lines == "" | str_detect(tab_lines, "^#"))
  if (length(end_idx) > 0) tab_lines <- tab_lines[seq_len(min(end_idx)-1)]

  # read into data.frame
  con <- textConnection(paste(tab_lines, collapse = "\n"))
  on.exit(close(con), add = TRUE)
  df <- read.table(con, header = FALSE, stringsAsFactors = FALSE)
  colnames(df) <- c("Codon","AA","Fraction","Frequency","Number")

  # drop stop codons (AA == "*")
  df <- df[df$AA != "*", ]
  df$Organism <- org_label
  df
}

# Read all organisms (no setwd; build absolute paths safely)
all_list <- vector("list", length(files))
for (i in seq_along(files)) {
  fp <- file.path(data_dir, files[i])
  if (!file.exists(fp)) stop("Missing file: ", fp)
  all_list[[i]] <- read_codon_table(fp, organism_labels[i])
}
codon_df <- do.call(rbind, all_list)

# Build the complete codon set (61 sense codons) and expand per organism
all_codons <- unique(codon_df[, c("Codon","AA")])
all_codons <- all_codons[order(all_codons$AA, all_codons$Codon), ]
orgs <- unique(codon_df$Organism)

expand_one <- function(org){
  d <- all_codons
  d$Organism <- org
  d <- d[, c("Organism","Codon","AA")]
  d
}
expanded <- do.call(rbind, lapply(orgs, expand_one))
expanded <- merge(expanded, codon_df[, c("Organism","Codon","AA","Fraction","Frequency","Number")],
                  by = c("Organism","Codon","AA"), all.x = TRUE)
# replace NAs with 0 in numeric columns
for (nm in c("Fraction","Frequency","Number")) {
  expanded[[nm]][is.na(expanded[[nm]])] <- 0
}

# Compute per-AA normalized frequency and RSCU (no dplyr needed)
aas <- unique(expanded$AA)
expanded$freq_within_AA <- 0
for (org in orgs) {
  for (aa in aas) {
    idx <- expanded$Organism == org & expanded$AA == aa
    s <- sum(expanded$Frequency[idx], na.rm = TRUE)
    if (s > 0) expanded$freq_within_AA[idx] <- expanded$Frequency[idx] / s
  }
}
expanded$RSCU <- 0
for (org in orgs) {
  for (aa in aas) {
    idx <- expanded$Organism == org & expanded$AA == aa
    total <- sum(expanded$Frequency[idx], na.rm = TRUE)
    n_syn <- sum(idx)
    if (n_syn > 0) {
      mean_exp <- total / n_syn
      if (mean_exp > 0) expanded$RSCU[idx] <- expanded$Frequency[idx] / mean_exp
    }
  }
}

# Build matrix: rows = organisms (pretty labels), cols = codons (grouped by AA)
codon_order <- all_codons$Codon
value_col <- if (use_rscu) "RSCU" else "freq_within_AA"

mat <- matrix(NA_real_,
              nrow = length(orgs),
              ncol = length(codon_order),
              dimnames = list(orgs, codon_order))

for (i in seq_along(orgs)) {
  org <- orgs[i]
  sub <- expanded[expanded$Organism == org, c("Codon", value_col)]
  sub <- sub[match(codon_order, sub$Codon), ]
  mat[i, ] <- sub[[value_col]]
}
# Replace rownames with your provided labels in original order
rownames(mat) <- organism_labels[match(rownames(mat), organism_labels)]

## ================== PLOTTING OPTIONS (VERTICAL MAIN HEATMAP) ==================
transpose_heatmap <- FALSE      # vertical: organisms on X, codons on Y
add_spacers       <- TRUE       # keep amino-acid group spacers
n_label_codons_for_full <- 25   # how many codon labels to show on the full heatmap
png_w <- 2600
png_h <- 1600
x_cex <- 0.9                    # organism label size (X)
y_cex <- 0.55                   # codon label size (Y)
x_rot <- 45                     # slight rotation for organism names
## ==============================================================================

# ----- insert NA spacer columns between amino-acid groups (optional) -----
mat_base <- mat  # keep a clean copy for variance/labels
if (add_spacers) {
  aa_vec <- all_codons$AA
  cods   <- all_codons$Codon
  blocks <- split(cods, aa_vec)
  with_spacers <- NULL
  for (i in seq_along(blocks)) {
    block <- blocks[[i]]
    block_mat <- mat[, block, drop = FALSE]
    with_spacers <- if (is.null(with_spacers)) block_mat else cbind(with_spacers, block_mat)
    if (i < length(blocks)) {
      sp <- matrix(NA_real_, nrow = nrow(mat), ncol = 1,
                   dimnames = list(rownames(mat), paste0("|", i)))
      with_spacers <- cbind(with_spacers, sp)
    }
  }
  mat <- with_spacers
}

# -------- Build the matrix to plot (vertical) --------
plot_mat <- if (transpose_heatmap) t(mat) else mat

# --- Pick N most-variable codons to LABEL (plot ALL codons) ---
var_cod <- apply(mat_base, 2, function(v) var(v, na.rm = TRUE))
label_codons <- names(sort(var_cod, decreasing = TRUE))[seq_len(min(n_label_codons_for_full, length(var_cod)))]

# --- Scales: organism labels on top, sparse codon labels on left ---
all_y <- rownames(plot_mat)                     # codons (+ spacers)
all_x <- colnames(plot_mat)                     # organisms (6)
y_is_sp <- substr(all_y, 1, 1) == "|"
y_keep  <- which(!y_is_sp & all_y %in% label_codons)
x_at    <- seq_len(ncol(plot_mat))              # label all organisms

scales_list <- list(
  x = list(
    at = x_at,
    labels = all_x,
    rot = 45,              # angled for readability
    cex = x_cex,
    alternating = 3        # show labels on the TOP instead of bottom
  ),
  y = list(
    at = y_keep,
    labels = all_y[y_keep],
    cex = y_cex
  )
)

# ---- Filenames ----
fname_base <- if (use_rscu) {
  sprintf("fig1_codon_heatmap_RSCU_vertical_topLabels_%d", n_label_codons_for_full)
} else {
  sprintf("fig1_codon_heatmap_freqWithinAA_vertical_topLabels_%d", n_label_codons_for_full)
}

# ---- PNG ----
png(file.path(out_dir, paste0(fname_base, ".png")), width = png_w, height = png_h, res = 140)
print(levelplot(
  plot_mat,
  xlab = "Organisms",
  ylab = "Codons (grouped by amino acid)",
  main = if (use_rscu)
    "RSCU heatmap (codon usage vs synonymous expectation)"
  else
    "Codon usage (frequency within amino acid)",
  scales = scales_list,
  col.regions = colorRampPalette(c("navy", "white", "firebrick"))(200),
  par.settings = list(fontsize = list(text = 9, points = 7)),
  xlim = extendrange(c(0.5, ncol(plot_mat) + 0.5))
))
dev.off()

# ---- PDF ----
pdf(file.path(out_dir, paste0(fname_base, ".pdf")), width = 14, height = 9)
print(levelplot(
  plot_mat,
  xlab = "Organisms",
  ylab = "Codons (grouped by amino acid)",
  main = if (use_rscu)
    "RSCU heatmap (codon usage vs synonymous expectation)"
  else
    "Codon usage (frequency within amino acid)",
  scales = scales_list,
  col.regions = colorRampPalette(c("navy", "white", "firebrick"))(200),
  par.settings = list(fontsize = list(text = 9, points = 7))
))
dev.off()

cat("✅ Full heatmap with top organism labels saved to:", normalizePath(out_dir), "\n")

# ----- (B) Optional: insert NA spacer columns between amino-acid groups -----
mat_base <- mat  # keep a clean copy for top-variance calc
if (add_spacers) {
  aa_vec <- all_codons$AA
  cods   <- all_codons$Codon
  blocks <- split(cods, aa_vec)

  with_spacers <- NULL
  for (i in seq_along(blocks)) {
    block <- blocks[[i]]
    block_mat <- mat[, block, drop = FALSE]
    if (is.null(with_spacers)) {
      with_spacers <- block_mat
    } else {
      with_spacers <- cbind(with_spacers, block_mat)
    }
    # add one NA spacer column after each AA block except the last
    if (i < length(blocks)) {
      sp <- matrix(NA_real_, nrow = nrow(mat), ncol = 1,
                   dimnames = list(rownames(mat), paste0("|", i)))
      with_spacers <- cbind(with_spacers, sp)
    }
  }
  mat <- with_spacers
}

# ----- (A) Plot the main heatmap (with transpose + bigger canvas) -----

## ---- OPTIONAL: Top-variance codons (keeps labels readable) ----
if (is.numeric(top_variance_N) && top_variance_N > 0) {
  cod_var <- apply(mat_base, 2, function(v) var(v, na.rm = TRUE))
  keep_codons <- names(sort(cod_var, decreasing = TRUE))[seq_len(min(top_variance_N, length(cod_var)))]
  mat_small <- mat_base[, keep_codons, drop = FALSE]

  # Thin X labels here too (codons along X after transpose below)
  pm <- t(mat_small)
  all_x <- colnames(pm)
  keepx <- seq(1, ncol(pm), by = every_n_codon_label)

  png(file.path(out_dir, sprintf("fig1b_topVariance_codons_%d.png", ncol(mat_small))),
      width = 1800, height = 1100, res = 140)
  print(levelplot(pm,
    xlab = "Organisms", ylab = sprintf("Top-%d variable codons", ncol(mat_small)),
    main = "Top-variance codons (compact view)",
    scales = list(
      x = list(cex = 0.95),
      y = list(cex = 0.75)  # codon labels on Y here
    ),
    col.regions = colorRampPalette(c("navy","white","firebrick"))(200),
    par.settings = list(fontsize = list(text = 9, points = 7))
  ))
  dev.off()
}