## ===== GO terms bubble plot (robust + factors) =====
library(lattice)

# ---- EDIT THESE ----
go_csv <- file.path("data", "go", "go_annotations_summary.csv")
out_dir <- file.path("figures", "generated_r")
# --------------------

stopifnot(file.exists(path.expand(go_csv)))
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# Read CSV; auto-skip a top title row if present
read_with_autoskip <- function(path){
  first <- readLines(path, n = 1, warn = FALSE)
  skip <- if (grepl(",", first, fixed = TRUE)) 0 else 1
  read.csv(path, stringsAsFactors = FALSE, check.names = FALSE,
           na.strings = c("", "NA"), skip = skip)
}
df_raw <- read_with_autoskip(go_csv)
stopifnot(nrow(df_raw) > 0)

# Clean headers/values
names(df_raw) <- trimws(names(df_raw))
for (j in seq_along(df_raw)) if (is.character(df_raw[[j]])) df_raw[[j]] <- trimws(df_raw[[j]])

need <- c("Organism (Full Name)", "GO Term", "Codon Count")
if (!all(need %in% names(df_raw))) {
  stop("CSV must have columns: ", paste(need, collapse = " | "),
       "\nFound: ", paste(names(df_raw), collapse = " | "))
}

# Tidy table
Count_num <- suppressWarnings(as.numeric(gsub("[^0-9.-]", "", as.character(df_raw[["Codon Count"]]))))
df_go <- data.frame(
  Organism = trimws(as.character(df_raw[["Organism (Full Name)"]])),
  GO       = trimws(as.character(df_raw[["GO Term"]])),
  Count    = Count_num,
  stringsAsFactors = FALSE
)
df_go <- df_go[nzchar(df_go$Organism) & nzchar(df_go$GO) & is.finite(df_go$Count), , drop = FALSE]

# Within-organism normalization + z-score
df_go$Percent <- ave(df_go$Count, df_go$Organism, FUN = function(x) {
  s <- sum(x, na.rm = TRUE); if (s > 0) x/s else NA_real_
})
df_go$Z <- ave(df_go$Percent, df_go$Organism, FUN = function(x) {
  mu <- mean(x, na.rm = TRUE); sdv <- sd(x, na.rm = TRUE)
  if (!is.finite(sdv) || sdv == 0) rep(0, length(x)) else (x - mu)/sdv
})

# Order GO labels by global prominence
ord <- sort(tapply(df_go$Percent, df_go$GO, median, na.rm = TRUE), decreasing = TRUE)
df_go$Label <- factor(df_go$GO, levels = rev(names(ord)))

# ---- FIX: make axes categorical factors (prevents NA coercion) ----
org_order <- c("Escherichia coli K-12",
               "Bacillus subtilis 168",
               "Salmonella enterica LT2",
               "Pseudomonas putida KT2440",
               "Halobacterium salinarum NRC-1",
               "Saccharomyces cerevisiae S288C")
present <- intersect(org_order, unique(df_go$Organism))
df_go$Organism <- factor(df_go$Organism,
                         levels = if (length(present)) org_order else unique(df_go$Organism))
df_go$Label <- factor(as.character(df_go$Label), levels = levels(df_go$Label))  # ensure non-NA factor

# Aesthetics: size ~ Percent, color ~ Z
z_to_col <- function(z) {
  pal <- colorRampPalette(c("navy","white","firebrick"))(100)
  zc  <- pmax(pmin(z, 3), -3)
  idx <- round((zc + 3) / 6 * 99) + 1
  pal[idx]
}
pt_cex <- pmax(0.8, 7 * sqrt(pmax(df_go$Percent, 0)))  # visible minimum
pt_col <- z_to_col(df_go$Z)

x_scales <- list(cex = 0.95, rot = 30)  # rotate organism labels a bit
y_scales <- list(cex = 0.7)

# ---- Bubble plot ----
png(file.path(out_dir, "fig2_GO_bubble.png"), width = 1700, height = 1200, res = 140)
print(xyplot(Label ~ Organism, data = df_go,
             pch = 16, cex = pt_cex, col = pt_col,
             xlab = "Organisms",
             ylab = "GO terms (union of top lists)",
             main = "GO prominence (size = % within organism; color = z-score within organism)",
             panel = function(x, y, subscripts, ...) {
               panel.grid(h = -1, v = -1, col = "grey85", lty = 3)
               panel.xyplot(x, y, subscripts = subscripts, ...)
             },
             scales = list(y = y_scales, x = x_scales)))
dev.off()

pdf(file.path(out_dir, "fig2_GO_bubble.pdf"), width = 12, height = 9)
print(xyplot(Label ~ Organism, data = df_go,
             pch = 16, cex = pt_cex, col = pt_col,
             xlab = "Organisms",
             ylab = "GO terms (union of top lists)",
             main = "GO prominence (size = % within organism; color = z-score within organism)",
             panel = function(x, y, subscripts, ...) {
               panel.grid(h = -1, v = -1, col = "grey85", lty = 3)
               panel.xyplot(x, y, subscripts = subscripts, ...)
             },
             scales = list(y = y_scales, x = x_scales)))
dev.off()

cat("✅ GO bubble plot saved to:", normalizePath(out_dir), "\n")

# ---- Top-N stacked bar per organism ----
topN <- 8
df_top <- do.call(rbind, by(df_go, df_go$Organism, function(dd) {
  o <- order(dd$Percent, decreasing = TRUE)
  dd[o[seq_len(min(topN, nrow(dd)))], ]
}))
labs <- levels(df_go$Label)
col_map <- setNames(colorRampPalette(c("#4575b4","#fee090","#d73027"))(length(labs)), labs)

png(file.path(out_dir, sprintf("fig2_GO_top%d_stacked.png", topN)), width = 1600, height = 900, res = 140)
print(barchart(Percent ~ Organism, groups = Label, data = df_top,
               stack = TRUE, horiz = FALSE, col = col_map[as.character(df_top$Label)],
               xlab = "Organisms", ylab = "Share within organism (sum of top-N)",
               main = sprintf("Top-%d GO terms per organism (stacked share)", topN),
               auto.key = list(columns = 3, cex = 0.8),
               scales = list(x = x_scales)))
dev.off()
