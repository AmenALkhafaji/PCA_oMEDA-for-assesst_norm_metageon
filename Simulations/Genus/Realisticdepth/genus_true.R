rm(list = ls())

# -------------------------
# Settings
# -------------------------
set.seed(1)

R <- 100
num_samples <- 150

# Depth ranges: realistic-depth balanced scenario
min_depth_healthy <- 10000
max_depth_healthy <- 100000

min_depth_unhealthy <- 10000
max_depth_unhealthy <- 100000

genera <- c(
  "Bifidobacterium", "Butyricimonas", "Odoribacter", "Paraprevotella",
  "Bacteroides", "Parabacteroides", "Prevotella", "Unknown1", "Unknown2",
  "Unknown3", "Clostridium", "Unknown4", "Unknown5", "Ruminococcus",
  "Anaerostipes", "Blautia", "Coprococcus", "Dorea", "Lachnospira",
  "Roseburia", "Unknown6", "Unknown7", "Unknown8", "Unknown9",
  "Faecalibacterium", "Oscillospira", "Ruminococcus_dup", "Dialister",
  "Unknown10", "Unknown11", "Unknown12", "Sutterella", "Escherichia",
  "Klebsiella", "Haemophilus", "Akkermansia"
)

prob_healthy <- c(
  0.15, 0.6, 1.5, 1.6, 12.2, 4.1, 30.8, 2.3, 0.14, 1.21,
  0.18, 0.32, 3.31, 0.63, 0.17, 0.34, 1.08, 0.14, 2.15,
  2.37, 0.033, 0.47, 0.34, 5.46, 3.47, 6.42, 4.25, 2.81,
  0.36, 1.02, 0.83, 1.8, 0.31, 1e-06, 0.17, 7.1
)

prob_unhealthy <- c(
  0.38, 0.6, 1.5, 0, 26, 0.01, 0.02, 11, 0.7, 1,
  0.14, 1, 0.15, 0.3, 0, 0.3, 3.00, 0.14, 0.02, 0.5,
  1, 0.14, 2, 4.8, 0, 0.2, 0.05, 5, 3, 2.7, 0.8,
  8, 0.22, 0, 0.31, 1.23
)


# -------------------------
# Helper: simulate one group
# -------------------------
simulate_group <- function(n, probs, group_label, min_depth, max_depth) {
  depths <- sample(min_depth:max_depth, size = n, replace = TRUE)
  
  mat <- t(sapply(depths, function(d) {
    rmultinom(n = 1, size = d, prob = probs)
  }))
  
  colnames(mat) <- genera
  
  df <- as.data.frame(mat)
  df$group <- factor(group_label)
  df
}

# -------------------------
# Output folders
# -------------------------
out_dir_counts <- "d:/Genus/True/RawData"
dir.create(out_dir_counts, showWarnings = FALSE, recursive = TRUE)

out_dir_rar <- "d:/Genus/True/Rarefaction"
dir.create(out_dir_rar, showWarnings = FALSE, recursive = TRUE)

out_dir_clr <- "d:/Genus/True/CLR"
dir.create(out_dir_clr, showWarnings = FALSE, recursive = TRUE)

out_dir_clr_mr <- "d:/Genus/True/CLR_MR"
dir.create(out_dir_clr_mr, showWarnings = FALSE, recursive = TRUE)

out_dir_clr_bmr <- "d:/Genus/True/CLR_BMR_safe"
dir.create(out_dir_clr_bmr, showWarnings = FALSE, recursive = TRUE)

out_dir_tss <- "d:/Genus/True/TSS"
dir.create(out_dir_tss, showWarnings = FALSE, recursive = TRUE)

out_dir_css <- "d:/Genus/True/CSS"
dir.create(out_dir_css, showWarnings = FALSE, recursive = TRUE)

out_dir_tmm <- "d:/Genus/True/edgeR_TMM"
dir.create(out_dir_tmm, showWarnings = FALSE, recursive = TRUE)

out_dir_deseq <- "d:/Genus/True/DESeq2"
dir.create(out_dir_deseq, showWarnings = FALSE, recursive = TRUE)


out_dir_aldex <- "d:/Genus/True/ALDEX"
dir.create(out_dir_aldex, showWarnings = FALSE, recursive = TRUE)

out_dir_ancom <- "d:/Genus/True/ANCOM"
dir.create(out_dir_ancom, showWarnings = FALSE, recursive = TRUE)

# -------------------------
# Rarefaction
# -------------------------
rarefy_only_taxa <- function(df, taxa_cols, depth = NULL) {
  if (!requireNamespace("vegan", quietly = TRUE)) {
    stop("Package 'vegan' is required. Install it with: install.packages('vegan')")
  }
  
  X <- as.matrix(df[, taxa_cols, drop = FALSE])
  
  if (!is.numeric(X)) stop("All taxa columns must be numeric.")
  if (any(is.na(X))) stop("Taxa matrix contains NA values.")
  if (any(X < 0)) stop("Counts must be non-negative.")
  
  libsize <- rowSums(X)
  
  if (any(libsize == 0)) {
    stop("At least one sample has total count zero; rarefaction cannot be computed.")
  }
  
  if (is.null(depth)) {
    depth <- min(libsize)
  }
  
  X_rar <- vegan::rrarefy(X, sample = depth)
  
  out <- as.data.frame(X_rar)
  colnames(out) <- taxa_cols
  rownames(out) <- rownames(df)
  out
}

# -------------------------
# CLR + pseudocount 1
# CLR+1
# -------------------------
clr_only_taxa <- function(df, taxa_cols, pseudocount = 1) {
  X <- as.matrix(df[, taxa_cols, drop = FALSE])
  
  if (!is.numeric(X)) stop("All taxa columns must be numeric.")
  if (any(is.na(X))) stop("Taxa matrix contains NA values.")
  if (any(X < 0)) stop("Counts must be non-negative.")
  
  X <- X + pseudocount
  
  logX <- log(X)
  clrX <- logX - rowMeans(logX)
  
  out <- as.data.frame(clrX)
  colnames(out) <- taxa_cols
  rownames(out) <- rownames(df)
  out
}

# -------------------------
# Helper: standard CLR on strictly positive proportions
# -------------------------
clr_from_positive_proportions <- function(X_prop, taxa_cols, row_ids = NULL) {
  X_prop <- as.matrix(X_prop)
  
  if (any(is.na(X_prop))) {
    stop("CLR input contains NA values.")
  }
  
  if (any(X_prop <= 0)) {
    stop("CLR on proportions requires strictly positive values.")
  }
  
  if (ncol(X_prop) != length(taxa_cols)) {
    stop(
      "Column mismatch before CLR: expected ",
      length(taxa_cols), " columns but got ", ncol(X_prop), "."
    )
  }
  
  logX <- log(X_prop)
  clrX <- logX - rowMeans(logX)
  
  out <- as.data.frame(clrX)
  colnames(out) <- taxa_cols
  
  if (!is.null(row_ids)) {
    rownames(out) <- row_ids
  }
  
  out
}

# -------------------------
# Helper: multiplicative replacement using zCompositions
# -------------------------
czm_replace <- function(X_prop) {
  X_repl <- zCompositions::cmultRepl(
    X_prop,
    label = 0,
    method = "CZM",
    output = "prop",
    z.delete = FALSE,
    z.warning = FALSE
  )
  
  X_repl <- as.matrix(X_repl)
  
  if (any(is.na(X_repl))) {
    stop("CZM replacement produced NA values.")
  }
  
  if (any(X_repl <= 0)) {
    stop("CZM replacement produced zero or negative values.")
  }
  
  X_repl <- sweep(X_repl, 1, rowSums(X_repl), "/")
  X_repl
}

# -------------------------
# CLR with multiplicative replacement
# CLR-MR
# -------------------------
clr_mr_only_taxa <- function(df, taxa_cols) {
  if (!requireNamespace("zCompositions", quietly = TRUE)) {
    stop("Package 'zCompositions' is required. Install it with: install.packages('zCompositions')")
  }
  
  X <- as.matrix(df[, taxa_cols, drop = FALSE])
  
  if (!is.numeric(X)) stop("All taxa columns must be numeric.")
  if (any(is.na(X))) stop("Taxa matrix contains NA values.")
  if (any(X < 0)) stop("Counts must be non-negative.")
  
  libsize <- rowSums(X)
  
  if (any(libsize == 0)) {
    stop("At least one sample has total count zero; CLR cannot be computed.")
  }
  
  X_prop <- sweep(X, 1, libsize, "/")
  
  if (!any(X_prop == 0)) {
    return(clr_from_positive_proportions(
      X_prop = X_prop,
      taxa_cols = taxa_cols,
      row_ids = rownames(df)
    ))
  }
  
  X_repl <- czm_replace(X_prop)
  
  if (ncol(X_repl) != length(taxa_cols)) {
    stop(
      "Column mismatch after multiplicative replacement: expected ",
      length(taxa_cols), " columns but got ", ncol(X_repl), "."
    )
  }
  
  clr_from_positive_proportions(
    X_prop = X_repl,
    taxa_cols = taxa_cols,
    row_ids = rownames(df)
  )
}

# -------------------------
# CLR with Bayesian-multiplicative replacement where possible
# CLR-BMR-safe
# -------------------------
clr_bmr_only_taxa <- function(df, taxa_cols) {
  if (!requireNamespace("zCompositions", quietly = TRUE)) {
    stop("Package 'zCompositions' is required. Install it with: install.packages('zCompositions')")
  }
  
  X <- as.matrix(df[, taxa_cols, drop = FALSE])
  
  if (!is.numeric(X)) stop("All taxa columns must be numeric.")
  if (any(is.na(X))) stop("Taxa matrix contains NA values.")
  if (any(X < 0)) stop("Counts must be non-negative.")
  
  libsize <- rowSums(X)
  
  if (any(libsize == 0)) {
    stop("At least one sample has total count zero; CLR cannot be computed.")
  }
  
  X_prop <- sweep(X, 1, libsize, "/")
  
  if (!any(X_prop == 0)) {
    return(clr_from_positive_proportions(
      X_prop = X_prop,
      taxa_cols = taxa_cols,
      row_ids = rownames(df)
    ))
  }
  
  positive_counts <- colSums(X_prop > 0)
  
  if (any(positive_counts < 2)) {
    message(
      "CLR-BMR-safe: at least one taxon has fewer than 2 positive values. ",
      "GBM cannot estimate its hyper-parameter reliably. ",
      "Using CZM fallback to preserve all taxa."
    )
    
    X_repl <- czm_replace(X_prop)
    
  } else {
    X_repl <- tryCatch(
      {
        zCompositions::cmultRepl(
          X_prop,
          label = 0,
          method = "GBM",
          output = "prop",
          z.delete = FALSE,
          z.warning = FALSE
        )
      },
      error = function(e) {
        message(
          "CLR-BMR-safe: GBM failed with message: ",
          conditionMessage(e),
          ". Using CZM fallback to preserve all taxa."
        )
        
        czm_replace(X_prop)
      }
    )
    
    X_repl <- as.matrix(X_repl)
    
    if (any(is.na(X_repl))) {
      stop("GBM/CZM replacement produced NA values.")
    }
    
    if (any(X_repl <= 0)) {
      stop("GBM/CZM replacement produced zero or negative values.")
    }
    
    X_repl <- sweep(X_repl, 1, rowSums(X_repl), "/")
  }
  
  if (ncol(X_repl) != length(taxa_cols)) {
    stop(
      "Column mismatch after CLR-BMR-safe replacement: expected ",
      length(taxa_cols), " columns but got ", ncol(X_repl), "."
    )
  }
  
  clr_from_positive_proportions(
    X_prop = X_repl,
    taxa_cols = taxa_cols,
    row_ids = rownames(df)
  )
}

# -------------------------
# TSS
# -------------------------
tss_only_taxa <- function(df, taxa_cols) {
  X <- as.matrix(df[, taxa_cols, drop = FALSE])
  
  if (!is.numeric(X)) stop("All taxa columns must be numeric.")
  if (any(is.na(X))) stop("Taxa matrix contains NA values.")
  if (any(X < 0)) stop("Counts must be non-negative.")
  
  libsize <- rowSums(X)
  libsize[libsize == 0] <- NA_real_
  
  X_tss <- X / libsize
  
  out <- as.data.frame(X_tss)
  colnames(out) <- taxa_cols
  rownames(out) <- rownames(df)
  out
}

# -------------------------
# CSS
# -------------------------
css_only_taxa <- function(df, taxa_cols) {
  if (!requireNamespace("metagenomeSeq", quietly = TRUE)) {
    stop("Package 'metagenomeSeq' is required. Install it with: BiocManager::install('metagenomeSeq')")
  }
  
  X <- as.matrix(df[, taxa_cols, drop = FALSE])
  
  if (!is.numeric(X)) stop("All taxa columns must be numeric.")
  if (any(is.na(X))) stop("Taxa matrix contains NA values.")
  if (any(X < 0)) stop("Counts must be non-negative.")
  
  X_t <- t(X)
  
  obj <- metagenomeSeq::newMRexperiment(X_t)
  p <- metagenomeSeq::cumNormStatFast(obj)
  obj <- metagenomeSeq::cumNorm(obj, p = p)
  
  X_css <- metagenomeSeq::MRcounts(obj, norm = TRUE, log = FALSE)
  
  out <- as.data.frame(t(X_css))
  colnames(out) <- taxa_cols
  rownames(out) <- rownames(df)
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
  
  if (!is.numeric(X)) stop("All taxa columns must be numeric.")
  if (any(is.na(X))) stop("Taxa matrix contains NA values.")
  if (any(X < 0)) stop("Counts must be non-negative.")
  
  dge <- edgeR::DGEList(counts = t(X))
  dge <- edgeR::calcNormFactors(dge, method = "TMM")
  tmm_cpm <- edgeR::cpm(dge, log = FALSE)
  
  out <- as.data.frame(t(tmm_cpm))
  colnames(out) <- taxa_cols
  rownames(out) <- rownames(df)
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
    stop(
      "Some taxa values became NA after numeric conversion. ",
      "This usually happens if the file was read incorrectly. ",
      "Use read.csv2 for csv2 files, or write.csv/read.csv for standard CSV."
    )
  }
  
  X <- as.matrix(Xdf)
  n_samp <- nrow(X)
  
  taxa_ids <- taxa_cols
  sample_ids <- paste0("S", seq_len(n_samp))
  
  countData <- t(round(X))
  rownames(countData) <- taxa_ids
  colnames(countData) <- sample_ids
  storage.mode(countData) <- "integer"
  
  colData <- data.frame(
    row.names = sample_ids,
    dummy = rep(1, n_samp)
  )
  
  stopifnot(ncol(countData) == nrow(colData))
  
  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = countData,
    colData = colData,
    design = ~ 1
  )
  
  dds <- DESeq2::estimateSizeFactors(dds, type = sf_type)
  norm_counts <- DESeq2::counts(dds, normalized = TRUE)
  
  out <- as.data.frame(t(norm_counts))
  colnames(out) <- taxa_cols
  rownames(out) <- rownames(df)
  out
}


# -------------------------
# ALDEx2-derived transformed output
# ALDEX = median CLR across ALDEx2 Monte Carlo instances
# This exports a transformed matrix for PCA/oMEDA; it is not a DA result.
#
# Important:
# ALDEx2 may internally drop taxa that are zero across all samples.
# This function returns the original taxa columns by inserting dropped taxa
# back as zero-valued CLR columns, so downstream matrices keep identical width.
# -------------------------
aldex_clr_median_only_taxa <- function(
    df,
    taxa_cols,
    group_col = "group",
    mc.samples = 128,
    denom = "all"
) {
  if (!requireNamespace("ALDEx2", quietly = TRUE)) {
    stop("Package 'ALDEx2' is required. Install it with: BiocManager::install('ALDEx2')")
  }
  
  X <- as.matrix(df[, taxa_cols, drop = FALSE])
  
  if (!is.numeric(X)) stop("All taxa columns must be numeric.")
  if (any(is.na(X))) stop("Taxa matrix contains NA values.")
  if (any(X < 0)) stop("Counts must be non-negative.")
  if (any(rowSums(X) == 0)) {
    stop("At least one sample has total count zero; ALDEx2 cannot be computed.")
  }
  
  if (!group_col %in% colnames(df)) {
    stop("Column '", group_col, "' not found in df.")
  }
  
  # ALDEx2 can fail or drop taxa that are zero across all samples.
  # We therefore run ALDEx2 only on taxa with at least one positive count.
  keep_taxa <- taxa_cols[colSums(X) > 0]
  dropped_taxa <- setdiff(taxa_cols, keep_taxa)
  
  if (length(keep_taxa) < 2) {
    stop("Fewer than two non-zero taxa are available for ALDEx2.")
  }
  
  if (length(dropped_taxa) > 0) {
    message(
      "ALDEX: dropping all-zero taxa before ALDEx2 and reinserting them as zero columns: ",
      paste(dropped_taxa, collapse = ", ")
    )
  }
  
  X_keep <- X[, keep_taxa, drop = FALSE]
  
  sample_ids <- paste0("S", seq_len(nrow(X_keep)))
  
  # ALDEx2 expects taxa/features in rows and samples in columns.
  X_t <- t(round(X_keep))
  rownames(X_t) <- keep_taxa
  colnames(X_t) <- sample_ids
  
  conds <- as.character(df[[group_col]])
  
  clr_obj <- tryCatch(
    {
      ALDEx2::aldex.clr(
        reads = X_t,
        conds = conds,
        mc.samples = mc.samples,
        denom = denom,
        verbose = FALSE,
        useMC = FALSE
      )
    },
    error = function(e) {
      stop("ALDEx2::aldex.clr failed: ", conditionMessage(e))
    }
  )
  
  if ("getMonteCarloInstances" %in% getNamespaceExports("ALDEx2")) {
    mc <- ALDEx2::getMonteCarloInstances(clr_obj)
  } else if (methods::is(clr_obj, "S4") && "analysisData" %in% slotNames(clr_obj)) {
    mc <- methods::slot(clr_obj, "analysisData")
  } else {
    stop("Could not extract ALDEx2 Monte Carlo instances.")
  }
  
  median_taxa_from_matrix <- function(m, keep_taxa) {
    m <- as.matrix(m)
    
    # Prefer names when available.
    if (!is.null(colnames(m))) {
      common <- intersect(keep_taxa, colnames(m))
      if (length(common) > 0) {
        vals <- rep(0, length(keep_taxa))
        names(vals) <- keep_taxa
        vals[common] <- apply(m[, common, drop = FALSE], 2, median, na.rm = TRUE)
        return(vals)
      }
    }
    
    if (!is.null(rownames(m))) {
      common <- intersect(keep_taxa, rownames(m))
      if (length(common) > 0) {
        vals <- rep(0, length(keep_taxa))
        names(vals) <- keep_taxa
        vals[common] <- apply(m[common, , drop = FALSE], 1, median, na.rm = TRUE)
        return(vals)
      }
    }
    
    # Fall back to dimensions.
    # Your observed failing case was 35 x 128:
    # 35 retained taxa x 128 Monte Carlo instances.
    if (nrow(m) == length(keep_taxa)) {
      vals <- apply(m, 1, median, na.rm = TRUE)
      names(vals) <- keep_taxa
      return(vals)
    }
    
    if (ncol(m) == length(keep_taxa)) {
      vals <- apply(m, 2, median, na.rm = TRUE)
      names(vals) <- keep_taxa
      return(vals)
    }
    
    stop(
      "Cannot infer taxa dimension in an ALDEx2 Monte Carlo matrix with dim: ",
      paste(dim(m), collapse = " x "),
      ". Expected one dimension to match retained taxa length = ",
      length(keep_taxa), "."
    )
  }
  
  if (is.list(mc) && length(mc) == nrow(X_keep)) {
    # Common structure: list by sample; each element is retained_taxa x MC or MC x retained_taxa.
    clr_med_keep <- sapply(mc, median_taxa_from_matrix, keep_taxa = keep_taxa)
    
  } else if (is.list(mc) && length(mc) == mc.samples) {
    # Alternative: list by Monte Carlo instance; each element is taxa x samples or samples x taxa.
    arr_list <- lapply(mc, as.matrix)
    first <- arr_list[[1]]
    
    if (nrow(first) == length(keep_taxa) && ncol(first) == nrow(X_keep)) {
      arr <- simplify2array(arr_list)
      clr_med_keep <- apply(arr, c(1, 2), median, na.rm = TRUE)
    } else if (nrow(first) == nrow(X_keep) && ncol(first) == length(keep_taxa)) {
      arr <- simplify2array(arr_list)
      clr_med_keep <- t(apply(arr, c(1, 2), median, na.rm = TRUE))
    } else {
      stop("Cannot infer ALDEx2 list-by-MC-instance dimensions.")
    }
    
  } else if (is.array(mc) && length(dim(mc)) == 3) {
    d <- dim(mc)
    
    if (d[1] == length(keep_taxa) && d[2] == nrow(X_keep)) {
      clr_med_keep <- apply(mc, c(1, 2), median, na.rm = TRUE)
    } else if (d[1] == nrow(X_keep) && d[2] == length(keep_taxa)) {
      clr_med_keep <- t(apply(mc, c(1, 2), median, na.rm = TRUE))
    } else if (d[2] == length(keep_taxa) && d[3] == nrow(X_keep)) {
      clr_med_keep <- apply(mc, c(2, 3), median, na.rm = TRUE)
    } else {
      stop("Unexpected ALDEx2 array dimensions: ", paste(d, collapse = " x "))
    }
    
  } else {
    stop(
      "Unexpected ALDEx2 Monte Carlo object structure. ",
      "Use str(ALDEx2::getMonteCarloInstances(clr_obj)) to inspect it."
    )
  }
  
  clr_med_keep <- as.matrix(clr_med_keep)
  
  if (nrow(clr_med_keep) != length(keep_taxa) || ncol(clr_med_keep) != nrow(X_keep)) {
    stop(
      "ALDEX retained-taxa output dimension mismatch: expected ",
      length(keep_taxa), " taxa x ", nrow(X_keep), " samples, got ",
      nrow(clr_med_keep), " x ", ncol(clr_med_keep), "."
    )
  }
  
  rownames(clr_med_keep) <- keep_taxa
  colnames(clr_med_keep) <- rownames(df)
  
  # Reinsert all original taxa columns.
  clr_med_full <- matrix(
    0,
    nrow = length(taxa_cols),
    ncol = nrow(X),
    dimnames = list(taxa_cols, rownames(df))
  )
  
  clr_med_full[keep_taxa, ] <- clr_med_keep[keep_taxa, , drop = FALSE]
  
  out <- as.data.frame(t(clr_med_full))
  colnames(out) <- taxa_cols
  rownames(out) <- rownames(df)
  
  out
}

# -------------------------
# ANCOM-BC-derived transformed output
# ANCOM = bias-corrected log observed abundance
# This exports a transformed matrix for PCA/oMEDA; it is not a DA result.
#
# Important:
# ANCOMBC has different ancombc() interfaces across package versions.
# This function first tries the newer "data = matrix/list" interface.
# If your installed version rejects taxa_are_rows/meta_data, it falls back
# to the older phyloseq interface.
# -------------------------
ancombc_logcorr_only_taxa <- function(
    df,
    taxa_cols,
    group_col = "group",
    pseudocount = 1
) {
  if (!requireNamespace("ANCOMBC", quietly = TRUE)) {
    stop("Package 'ANCOMBC' is required. Install it with: BiocManager::install('ANCOMBC')")
  }
  
  X <- as.matrix(df[, taxa_cols, drop = FALSE])
  
  if (!is.numeric(X)) stop("All taxa columns must be numeric.")
  if (any(is.na(X))) stop("Taxa matrix contains NA values.")
  if (any(X < 0)) stop("Counts must be non-negative.")
  if (any(rowSums(X) == 0)) {
    stop("At least one sample has total count zero; ANCOM-BC cannot be computed.")
  }
  
  if (!group_col %in% colnames(df)) {
    stop("Column '", group_col, "' not found in df.")
  }
  
  sample_ids <- paste0("S", seq_len(nrow(X)))
  
  # ANCOMBC generally expects taxa/features in rows and samples in columns.
  X_t <- t(round(X))
  rownames(X_t) <- taxa_cols
  colnames(X_t) <- sample_ids
  
  meta_data <- data.frame(
    group = factor(df[[group_col]]),
    row.names = sample_ids
  )
  
  out_ancom <- NULL
  
  # Try interface used by newer ANCOMBC-related functions.
  out_ancom <- tryCatch(
    {
      ANCOMBC::ancombc(
        data = X_t,
        taxa_are_rows = TRUE,
        meta_data = meta_data,
        formula = "group",
        p_adj_method = "BH",
        prv_cut = 0,
        lib_cut = 0,
        group = NULL,
        struc_zero = FALSE,
        neg_lb = FALSE,
        tol = 1e-5,
        max_iter = 100,
        conserve = TRUE,
        alpha = 0.05,
        global = FALSE,
        n_cl = 1,
        verbose = FALSE
      )
    },
    error = function(e1) {
      msg1 <- conditionMessage(e1)
      
      # Fall back for older ANCOMBC versions:
      # older ancombc() often expects a phyloseq object via phyloseq = ...
      if (!grepl("unused argument", msg1, fixed = TRUE) &&
          !grepl("unused arguments", msg1, fixed = TRUE)) {
        stop("ANCOMBC::ancombc failed: ", msg1)
      }
      
      if (!requireNamespace("phyloseq", quietly = TRUE)) {
        stop(
          "ANCOMBC::ancombc rejected the newer arguments, and package 'phyloseq' ",
          "is required for the older ANCOMBC interface. Install it with: ",
          "BiocManager::install('phyloseq'). Original ANCOMBC error: ", msg1
        )
      }
      
      message(
        "ANCOM: installed ANCOMBC uses an older interface. ",
        "Retrying with phyloseq = physeq."
      )
      
      physeq <- phyloseq::phyloseq(
        phyloseq::otu_table(X_t, taxa_are_rows = TRUE),
        phyloseq::sample_data(meta_data)
      )
      
      tryCatch(
        {
          ANCOMBC::ancombc(
            phyloseq = physeq,
            formula = "group",
            p_adj_method = "BH",
            prv_cut = 0,
            lib_cut = 0,
            group = NULL,
            struc_zero = FALSE,
            neg_lb = FALSE,
            tol = 1e-5,
            max_iter = 100,
            conserve = TRUE,
            alpha = 0.05,
            global = FALSE,
            n_cl = 1,
            verbose = FALSE
          )
        },
        error = function(e2) {
          stop(
            "ANCOMBC::ancombc failed with both interfaces. ",
            "New-interface error: ", msg1,
            " | phyloseq-interface error: ", conditionMessage(e2)
          )
        }
      )
    }
  )
  
  # Extract sample-specific sampling fractions.
  if (!is.null(out_ancom$samp_frac)) {
    samp_frac <- out_ancom$samp_frac
  } else if (!is.null(out_ancom$bias_correct_log_table) &&
             !is.null(out_ancom$feature_table)) {
    # Some versions may not expose samp_frac cleanly.
    # In that case, use the bias-corrected log table directly below.
    samp_frac <- NULL
  } else {
    stop(
      "ANCOM-BC output does not contain samp_frac. Available output names are: ",
      paste(names(out_ancom), collapse = ", ")
    )
  }
  
  # Preferred direct output if available.
  if (!is.null(out_ancom$bias_correct_log_table)) {
    log_corr_abn <- as.matrix(out_ancom$bias_correct_log_table)
    
    # Ensure taxa x samples orientation.
    if (all(taxa_cols %in% rownames(log_corr_abn))) {
      log_corr_abn <- log_corr_abn[taxa_cols, , drop = FALSE]
    } else if (all(taxa_cols %in% colnames(log_corr_abn))) {
      log_corr_abn <- t(log_corr_abn[, taxa_cols, drop = FALSE])
    } else {
      stop("ANCOM-BC bias_correct_log_table does not contain the expected taxa names.")
    }
    
    # Ensure sample order.
    if (all(sample_ids %in% colnames(log_corr_abn))) {
      log_corr_abn <- log_corr_abn[, sample_ids, drop = FALSE]
    }
    
  } else {
    # Reconstruct corrected log abundance:
    # log observed abundance minus estimated sample-specific sampling fraction.
    samp_frac[is.na(samp_frac)] <- 0
    
    if (!is.null(names(samp_frac))) {
      samp_frac <- samp_frac[sample_ids]
    }
    
    feature_table <- NULL
    
    if (!is.null(out_ancom$feature_table)) {
      feature_table <- as.matrix(out_ancom$feature_table)
    } else {
      feature_table <- X_t
    }
    
    # Ensure taxa x samples orientation.
    if (all(taxa_cols %in% rownames(feature_table))) {
      feature_table <- feature_table[taxa_cols, , drop = FALSE]
    } else if (all(taxa_cols %in% colnames(feature_table))) {
      feature_table <- t(feature_table[, taxa_cols, drop = FALSE])
    } else {
      stop("Some taxa are missing from ANCOM-BC feature_table.")
    }
    
    if (all(sample_ids %in% colnames(feature_table))) {
      feature_table <- feature_table[, sample_ids, drop = FALSE]
    }
    
    log_obs_abn <- log(feature_table + pseudocount)
    log_corr_abn <- t(t(log_obs_abn) - samp_frac)
  }
  
  # Final validation and export as samples x taxa.
  if (nrow(log_corr_abn) != length(taxa_cols)) {
    stop(
      "ANCOM output taxa dimension mismatch: expected ",
      length(taxa_cols), " taxa but got ", nrow(log_corr_abn), "."
    )
  }
  
  out <- as.data.frame(t(log_corr_abn))
  colnames(out) <- taxa_cols
  rownames(out) <- rownames(df)
  
  out
}

# -------------------------
# Helper: stratify each group by depth and interleave bands
# -------------------------
make_banded_sample_summary <- function(healthy_df, unhealthy_df, taxa_cols, block_size = 50) {
  healthy_df$depth_sum <- rowSums(healthy_df[, taxa_cols, drop = FALSE])
  unhealthy_df$depth_sum <- rowSums(unhealthy_df[, taxa_cols, drop = FALSE])
  
  healthy_sorted <- healthy_df[order(healthy_df$depth_sum, decreasing = FALSE), ]
  unhealthy_sorted <- unhealthy_df[order(unhealthy_df$depth_sum, decreasing = FALSE), ]
  
  if (nrow(healthy_sorted) < 3 * block_size || nrow(unhealthy_sorted) < 3 * block_size) {
    stop("Not enough samples to form low/mid/high bands of size ", block_size, " per group.")
  }
  
  H_low <- healthy_sorted[1:block_size, ]
  H_mid <- healthy_sorted[(block_size + 1):(2 * block_size), ]
  H_high <- healthy_sorted[(2 * block_size + 1):(3 * block_size), ]
  
  U_low <- unhealthy_sorted[1:block_size, ]
  U_mid <- unhealthy_sorted[(block_size + 1):(2 * block_size), ]
  U_high <- unhealthy_sorted[(2 * block_size + 1):(3 * block_size), ]
  
  sample_summary <- rbind(H_low, U_low, H_mid, U_mid, H_high, U_high)
  
  sample_summary$depth_sum <- NULL
  rownames(sample_summary) <- NULL
  
  sample_summary
}

# -------------------------
# Run repetitions and save
# -------------------------
for (r in 1:R) {
  
  healthy_df <- simulate_group(
    n = num_samples,
    probs = prob_healthy,
    group_label = 1,
    min_depth = min_depth_healthy,
    max_depth = max_depth_healthy
  )
  
  unhealthy_df <- simulate_group(
    n = num_samples,
    probs = prob_unhealthy,
    group_label = 2,
    min_depth = min_depth_unhealthy,
    max_depth = max_depth_unhealthy
  )
  
  sample_summary <- make_banded_sample_summary(
    healthy_df,
    unhealthy_df,
    genera,
    block_size = 50
  )
  
  zero_counts <- sum(sample_summary[, genera] == 0)
  zero_percentage <- zero_counts / (length(genera) * nrow(sample_summary)) * 100
  
  message(sprintf("Rep %03d: zero%% = %.3f", r, zero_percentage))
  
  # Raw counts
  file_counts <- file.path(
    out_dir_counts,
    sprintf("genus_counts_rep%03d.csv", r)
  )
  write.csv(sample_summary, file_counts, row.names = FALSE)
  
  # edgeR TMM-CPM
  tmm_5 <- edger_tmm_cpm_only_taxa(sample_summary, genera)
  file_tmm <- file.path(
    out_dir_tmm,
    sprintf("genus_edgeR_TMM_CPM_rep%03d.csv", r)
  )
  write.csv(tmm_5, file_tmm, row.names = FALSE)
  
  # DESeq2 normalized counts
  deseq_5 <- deseq_norm_only_taxa(sample_summary, genera, sf_type = "poscounts")
  file_deseq <- file.path(
    out_dir_deseq,
    sprintf("genus_DESeq2norm_rep%03d.csv", r)
  )
  write.csv(deseq_5, file_deseq, row.names = FALSE)
  
  # CSS
  css_5 <- css_only_taxa(sample_summary, genera)
  file_css <- file.path(
    out_dir_css,
    sprintf("genus_CSS_rep%03d.csv", r)
  )
  write.csv(css_5, file_css, row.names = FALSE)
  
  # TSS
  tss_5 <- tss_only_taxa(sample_summary, genera)
  file_tss <- file.path(
    out_dir_tss,
    sprintf("genus_TSS_rep%03d.csv", r)
  )
  write.csv(tss_5, file_tss, row.names = FALSE)
  
  # CLR+1: CLR with pseudocount 1
  clr_5 <- clr_only_taxa(sample_summary, genera, pseudocount = 1)
  file_clr <- file.path(
    out_dir_clr,
    sprintf("genus_CLR_pseudocount1_rep%03d.csv", r)
  )
  write.csv(clr_5, file_clr, row.names = FALSE)
  
  # CLR-MR: CLR after multiplicative replacement
  clr_mr_5 <- clr_mr_only_taxa(sample_summary, genera)
  file_clr_mr <- file.path(
    out_dir_clr_mr,
    sprintf("genus_CLR_MR_rep%03d.csv", r)
  )
  write.csv(clr_mr_5, file_clr_mr, row.names = FALSE)
  
  # CLR-BMR-safe: GBM when estimable; CZM fallback otherwise
  clr_bmr_5 <- clr_bmr_only_taxa(sample_summary, genera)
  file_clr_bmr <- file.path(
    out_dir_clr_bmr,
    sprintf("genus_CLR_BMR_safe_rep%03d.csv", r)
  )
  write.csv(clr_bmr_5, file_clr_bmr, row.names = FALSE)
  
  
  # ALDEX: ALDEx2-derived median CLR output
  aldex_5 <- aldex_clr_median_only_taxa(
    sample_summary,
    genera,
    group_col = "group",
    mc.samples = 128,
    denom = "all"
  )
  file_aldex <- file.path(
    out_dir_aldex,
    sprintf("genus_ALDEX_rep%03d.csv", r)
  )
  write.csv(aldex_5, file_aldex, row.names = FALSE)
  
  # ANCOM: ANCOM-BC-derived bias-corrected log-abundance output
  ancom_5 <- ancombc_logcorr_only_taxa(
    sample_summary,
    genera,
    group_col = "group",
    pseudocount = 1
  )
  file_ancom <- file.path(
    out_dir_ancom,
    sprintf("genus_ANCOM_rep%03d.csv", r)
  )
  write.csv(ancom_5, file_ancom, row.names = FALSE)
  
  # Rarefaction
  rarefied_5 <- rarefy_only_taxa(sample_summary, genera)
  file_rar <- file.path(
    out_dir_rar,
    sprintf("genus_RAREFY_rep%03d.csv", r)
  )
  write.csv2(rarefied_5, file_rar, row.names = FALSE)
}