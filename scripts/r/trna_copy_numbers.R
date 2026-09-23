## ================== tRNA copy-number scatter (from XLSX) ==================
library(lattice)

# --------- USER: set your XLSX path ----------
trna_csv <- file.path("data", "trnascan", "summary_outputs.csv")  
stopifnot(file.exists(trna_csv))

# Were to save
out_dir <- file.path("figures", "generated_r")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# --------- 1) Read & clean ----------
# Expect columns: Organism, Codon, Amino Acid, Count
df0 <- read.csv(trna_csv, check.names = FALSE, stringsAsFactors = FALSE)
nm  <- tolower(names(df0))

# try to auto-match headers that may vary slightly
col_org   <- which(grepl("organism", nm))
col_codon <- which(grepl("^codon$", nm) | grepl("anti?codon", nm))
col_aa    <- which(grepl("^amino", nm) | grepl("^aa$", nm))
col_cnt   <- which(grepl("count|copies?", nm))

if (any(vapply(list(col_org, col_codon, col_aa, col_cnt), length, integer(1)) != 1L)) {
  stop("Couldn't find expected columns. Headers found:\n",
       paste(names(df0), collapse = " | "),
       "\nNeed something like: Organism, Codon, Amino Acid, Count")
}

df <- data.frame(
  Organism  = as.character(df0[[col_org]]),
  Codon     = toupper(gsub("\\s+", "", as.character(df0[[col_codon]]))), # remove spaces in codons like 'G T T'
  AA        = as.character(df0[[col_aa]]),
  Count     = suppressWarnings(as.numeric(df0[[col_cnt]])),
  stringsAsFactors = FALSE
)

# drop empty/invalid
df <- df[ nzchar(df$Organism) & nzchar(df$Codon) & !is.na(df$Count), , drop = FALSE]

# some tRNAs are listed as Ile2, fMet etc — keep labels as provided
df$AA[df$AA == ""] <- "Unknown"

# Aggregate in case the sheet has duplicates per Organism/Codon/AA
agg <- aggregate(Count ~ Organism + Codon + AA, data = df, sum, na.rm = TRUE)

# --------- 2) Build codon order grouped by amino acid (with spacers) ----------
# Order amino acids alphabetically; within each AA, codons alphabetically
aa_levels   <- sort(unique(agg$AA))
blocks      <- lapply(aa_levels, function(a) sort(unique(agg$Codon[agg$AA == a])))
names(blocks) <- aa_levels

# Create a numeric x-position for each codon and add thin NA spacer between AA blocks
x_pos <- numeric(0)
x_lab <- character(0)
x_aa  <- character(0)
pos <- 0
for (aa in aa_levels) {
  cods <- blocks[[aa]]
  idxs <- seq_len(length(cods)) + pos
  x_pos <- c(x_pos, idxs)
  x_lab <- c(x_lab, cods)
  x_aa  <- c(x_aa,  rep(aa, length(cods)))
  pos <- max(idxs) + 1  # leave 1 slot as spacer
}
# map codon -> numeric position
codon_to_x <- setNames(x_pos, x_lab)

agg$x <- codon_to_x[agg$Codon]
agg$AA <- factor(agg$AA, levels = aa_levels)

# --------- 3) Plot settings ----------
# label thinning so x-axis doesn’t crowd
every_n <- 3
tick_at  <- x_pos[seq(1, length(x_pos), by = every_n)]
tick_lab <- x_lab[seq(1, length(x_lab), by = every_n)]

# color palette for amino acids (simple base palette)
pal <- grDevices::rainbow(length(aa_levels))
aa_col <- setNames(pal, aa_levels)

# Arrange facets (panels) in a 3x2 grid
n_org <- length(unique(agg$Organism))
lay <- c(3, 2)  # 3 columns x 2 rows

# --------- 4) Make scatter/dot plot ----------
png(file.path(out_dir, "fig3_tRNA_scatter_by_codon.png"), width = 1900, height = 1100, res = 140)
print(xyplot(Count ~ x | Organism, data = agg,
             groups = AA,
             type = "p",
             pch = 16,
             col = aa_col[agg$AA],
             cex = pmin(1.2, 0.6 + 0.10 * sqrt(agg$Count)),  # a touch bigger for higher counts
             layout = lay,
             xlab = "Codons (grouped by amino acid; thin gaps = AA boundaries)",
             ylab = "tRNA copy number (tRNAscan-SE)",
             main = "tRNA copy numbers per codon (one panel per organism)",
             panel = function(x, y, subscripts, groups, ...) {
               # light vertical guides
               panel.abline(v = tick_at, col = "grey90", lty = 3)
               panel.xyplot(x, y, subscripts = subscripts, ...)
             },
             auto.key = list(points = TRUE, columns = 5, title = "Amino acid", cex.title = 0.9),
             scales = list(
               x = list(at = tick_at, labels = tick_lab, cex = 0.7, rot = 90),
               y = list(cex = 0.9)
             )))
dev.off()

pdf(file.path(out_dir, "fig3_tRNA_scatter_by_codon.pdf"), width = 14, height = 8)
print(xyplot(Count ~ x | Organism, data = agg,
             groups = AA,
             type = "p",
             pch = 16,
             col = aa_col[agg$AA],
             cex = pmin(1.2, 0.6 + 0.10 * sqrt(agg$Count)),
             layout = lay,
             xlab = "Codons (grouped by amino acid; thin gaps = AA boundaries)",
             ylab = "tRNA copy number (tRNAscan-SE)",
             main = "tRNA copy numbers per codon (one panel per organism)",
             panel = function(x, y, subscripts, groups, ...) {
               panel.abline(v = tick_at, col = "grey90", lty = 3)
               panel.xyplot(x, y, subscripts = subscripts, ...)
             },
             auto.key = list(points = TRUE, columns = 5, title = "Amino acid", cex.title = 0.9),
             scales = list(
               x = list(at = tick_at, labels = tick_lab, cex = 0.7, rot = 90),
               y = list(cex = 0.9)
             )))
dev.off()

cat("✅ tRNA scatter saved to:", normalizePath(out_dir), "\n")
