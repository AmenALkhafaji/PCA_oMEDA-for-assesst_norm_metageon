rm(list = ls())

# -------------------------
# Settings
# -------------------------
set.seed(1)
R <- 100
num_samples <- 150
min_depth <- 100
max_depth <- 300

phyla <- c("Bacillota", "Bacteroidetes", "Actinobacteria", "Verrucomicrobiota", "Proteobacteria")

prob_healthy   <- c(50, 35, 0.5, 6, 2.2)

prob_unhealthy <- c(50, 35, 0.5, 6, 2.2)

# -------------------------
# Helper: simulate one group
# -------------------------
simulate_group <- function(n, probs, group_label) {
  depths <- sample(min_depth:max_depth, size = n, replace = TRUE)
  mat <- t(sapply(depths, function(d) rmultinom(n = 1, size = d, prob = probs)))
  colnames(mat) <- phyla
  df <- as.data.frame(mat)
  df$group <- factor(group_label)
  df
}

# -------------------------
# Helper: stratify each group by depth and interleave bands
# -------------------------
make_banded_sample_summary <- function(healthy_df, unhealthy_df, taxa_cols, block_size = 50) {
  healthy_df$depth_sum   <- rowSums(healthy_df[, taxa_cols, drop = FALSE])
  unhealthy_df$depth_sum <- rowSums(unhealthy_df[, taxa_cols, drop = FALSE])
  
  healthy_sorted   <- healthy_df[order(healthy_df$depth_sum, decreasing = FALSE), ]
  unhealthy_sorted <- unhealthy_df[order(unhealthy_df$depth_sum, decreasing = FALSE), ]
  
  if (nrow(healthy_sorted) < 3 * block_size || nrow(unhealthy_sorted) < 3 * block_size) {
    stop("Not enough samples to form low/mid/high bands of size ", block_size, " per group.")
  }
  
  H_low  <- healthy_sorted[1:block_size, ]
  H_mid  <- healthy_sorted[(block_size + 1):(2 * block_size), ]
  H_high <- healthy_sorted[(2 * block_size + 1):(3 * block_size), ]
  
  U_low  <- unhealthy_sorted[1:block_size, ]
  U_mid  <- unhealthy_sorted[(block_size + 1):(2 * block_size), ]
  U_high <- unhealthy_sorted[(2 * block_size + 1):(3 * block_size), ]
  
  sample_summary <- rbind(H_low, U_low, H_mid, U_mid, H_high, U_high)
  
  sample_summary$depth_sum <- NULL
  rownames(sample_summary) <- NULL
  sample_summary
}

# -------------------------
# Helper: TSS normalize
# -------------------------
tss_normalize <- function(df, taxa_cols) {
  X <- as.matrix(df[, taxa_cols, drop = FALSE])
  rs <- rowSums(X)
  rs[rs == 0] <- NA_real_
  
  X_tss <- X / rs
  out <- as.data.frame(X_tss)
  out$group <- df$group
  out
}

# -------------------------
# Output folders
# -------------------------
out_dir_counts <- "d:/True_sim_phylum_reps_counts"
out_dir_tss    <- "d:/True_sim_phylum_reps_TSS"
out_dir_rar    <- "d:/True_sim_phylum_reps_RAREFY"
out_dir_clr    <- "d:/True_sim_phylum_reps_CLR"
out_dir_css    <- "d:/True_sim_phylum_reps_CSS"
out_dir_tmm    <- "d:/True_sim_phylum_reps_edgeR_TMM_CPM"
out_dir_deseq  <- "d:/True_sim_phylum_reps_DESeq2"

dir.create(out_dir_counts, showWarnings = FALSE, recursive = TRUE)
dir.create(out_dir_tss,    showWarnings = FALSE, recursive = TRUE)
dir.create(out_dir_rar,    showWarnings = FALSE, recursive = TRUE)
dir.create(out_dir_clr,    showWarnings = FALSE, recursive = TRUE)
dir.create(out_dir_css,    showWarnings = FALSE, recursive = TRUE)
dir.create(out_dir_tmm,    showWarnings = FALSE, recursive = TRUE)
dir.create(out_dir_deseq,  showWarnings = FALSE, recursive = TRUE)

# -------------------------
# Rarefaction
# -------------------------
rarefy_only_taxa <- function(df, taxa_cols, depth = NULL) {
  if (!requireNamespace("vegan", quietly = TRUE)) {
    stop("Package 'vegan' is required. Install it with: install.packages('vegan')")
  }
  
  X <- as.matrix(df[, taxa_cols, drop = FALSE])
  libsize <- rowSums(X)
  
  if (is.null(depth)) depth <- min(libsize)
  
  X_rar <- vegan::rrarefy(X, sample = depth)
  as.data.frame(X_rar)
}

# -------------------------
# CLR
# -------------------------
clr_only_taxa <- function(df, taxa_cols, pseudocount = 1) {
  X <- as.matrix(df[, taxa_cols, drop = FALSE])
  X <- X + pseudocount
  logX <- log(X)
  clrX <- logX - rowMeans(logX)
  as.data.frame(clrX)
}

# -------------------------
# TSS
# -------------------------
tss_only_taxa <- function(df, taxa_cols) {
  X <- as.matrix(df[, taxa_cols, drop = FALSE])
  libsize <- rowSums(X)
  libsize[libsize == 0] <- NA_real_
  X_tss <- X / libsize
  as.data.frame(X_tss)
}

# -------------------------
# CSS
# -------------------------
css_only_taxa <- function(df, taxa_cols) {
  if (!requireNamespace("metagenomeSeq", quietly = TRUE)) {
    stop("Package 'metagenomeSeq' is required. Install it with: BiocManager::install('metagenomeSeq')")
  }
  
  X <- as.matrix(df[, taxa_cols, drop = FALSE])
  X_t <- t(X)
  
  obj <- metagenomeSeq::newMRexperiment(X_t)
  p <- metagenomeSeq::cumNormStatFast(obj)
  obj <- metagenomeSeq::cumNorm(obj, p = p)
  
  X_css <- metagenomeSeq::MRcounts(obj, norm = TRUE, log = FALSE)
  out <- as.data.frame(t(X_css))
  colnames(out) <- taxa_cols
  out
}

# -------------------------
# edgeR TMM-CPM
# -------------------------
edger_tmm_cpm_only_taxa <- function(df, taxa_cols) {
  if (!requireNamespace("edgeR", quietly = TRUE)) {
    stop("Package 'edgeR' is required. Install it with: BiocManager::install('edgeR')")
  }
  
  X <- as.matrix(df[, taxa_cols, drop = FALSE])
  
  if (!is.numeric(X)) {
    stop("All taxa columns must be numeric.")
  }
  
  if (any(is.na(X))) {
    stop("Taxa matrix contains NA values.")
  }
  
  if (any(X < 0)) {
    stop("Counts must be non-negative.")
  }
  
  dge <- edgeR::DGEList(counts = t(X))
  dge <- edgeR::calcNormFactors(dge, method = "TMM")
  tmm_cpm <- edgeR::cpm(dge, log = FALSE)
  
  out <- as.data.frame(t(tmm_cpm))
  colnames(out) <- taxa_cols
  out
}

# -------------------------
# DESeq2
# -------------------------
deseq_norm_only_taxa <- function(df, taxa_cols, sf_type = "poscounts") {
  if (!requireNamespace("DESeq2", quietly = TRUE)) {
    stop("Package 'DESeq2' is required. Install it with: BiocManager::install('DESeq2')")
  }
  
  Xdf <- df[, taxa_cols, drop = FALSE]
  
  for (j in seq_along(taxa_cols)) {
    Xdf[[j]] <- suppressWarnings(as.numeric(Xdf[[j]]))
  }
  
  if (anyNA(as.matrix(Xdf))) {
    stop("Some taxa values became NA after numeric conversion. ",
         "This usually happens if the file was read incorrectly (e.g., csv2/Excel separator). ",
         "Use read.csv2 for csv2 files, or write.csv/read.csv for standard CSV.")
  }
  
  X <- as.matrix(Xdf)
  
  n_samp <- nrow(X)
  n_tax  <- ncol(X)
  
  sample_ids <- paste0("S", seq_len(n_samp))
  taxa_ids   <- taxa_cols
  
  countData <- t(round(X))
  rownames(countData) <- taxa_ids
  colnames(countData) <- sample_ids
  storage.mode(countData) <- "integer"
  
  colData <- data.frame(row.names = sample_ids, dummy = rep(1, n_samp))
  stopifnot(ncol(countData) == nrow(colData))
  
  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = countData,
    colData   = colData,
    design    = ~ 1
  )
  
  dds <- DESeq2::estimateSizeFactors(dds, type = sf_type)
  norm_counts <- DESeq2::counts(dds, normalized = TRUE)
  
  out <- as.data.frame(t(norm_counts))
  colnames(out) <- taxa_cols
  out
}

# -------------------------
# Run repetitions and save
# -------------------------
for (r in 1:R) {
  healthy_df   <- simulate_group(num_samples, prob_healthy,   1)
  unhealthy_df <- simulate_group(num_samples, prob_unhealthy, 2)
  
  sample_summary <- make_banded_sample_summary(healthy_df, unhealthy_df, phyla, block_size = 50)
  
  zero_counts <- sum(sample_summary[, phyla] == 0)
  zero_percentage <- zero_counts / (length(phyla) * nrow(sample_summary)) * 100
  message(sprintf("Rep %03d: zero%% = %.3f", r, zero_percentage))
  
  # Raw counts
  file_counts <- file.path(out_dir_counts, sprintf("phylum_counts_rep%03d.csv", r))
  write.csv2(sample_summary, file_counts, row.names = FALSE)
  
  # Rarefaction
  rarefied_5 <- rarefy_only_taxa(sample_summary, phyla)
  file_rar <- file.path(out_dir_rar, sprintf("phylum_RAREFY_rep%03d.csv", r))
  write.csv2(rarefied_5, file_rar, row.names = FALSE)
  
  # CLR
  clr_5 <- clr_only_taxa(sample_summary, phyla, pseudocount = 1)
  file_clr <- file.path(out_dir_clr, sprintf("phylum_CLR_rep%03d.csv", r))
  write.csv(clr_5, file_clr, row.names = FALSE)
  
  # TSS
  tss_5 <- tss_only_taxa(sample_summary, phyla)
  file_tss <- file.path(out_dir_tss, sprintf("phylum_TSS_rep%03d.csv", r))
  write.csv(tss_5, file_tss, row.names = FALSE)
  
  # DESeq2
  deseq_5 <- deseq_norm_only_taxa(sample_summary, phyla, sf_type = "poscounts")
  file_deseq <- file.path(out_dir_deseq, sprintf("phylum_DESeq2norm_rep%03d.csv", r))
  write.csv(deseq_5, file_deseq, row.names = FALSE)
  
  # CSS
  css_5 <- css_only_taxa(sample_summary, phyla)
  file_css <- file.path(out_dir_css, sprintf("phylum_CSS_rep%03d.csv", r))
  write.csv(css_5, file_css, row.names = FALSE)
  
  # edgeR TMM-CPM
  tmm_5 <- edger_tmm_cpm_only_taxa(sample_summary, phyla)
  file_tmm <- file.path(out_dir_tmm, sprintf("phylum_edgeR_TMM_CPM_rep%03d.csv", r))
  write.csv(tmm_5, file_tmm, row.names = FALSE)
}