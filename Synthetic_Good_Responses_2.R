library(lavaan)
library(dplyr)

# ------------------------------------------------------------
# Helper: generate a valid factor correlation matrix
# Pairwise correlations are sampled from [min_cor, max_cor]
# and accepted only if the full matrix is positive definite.
# ------------------------------------------------------------
random_factor_cor_matrix <- function(nConstructs,
                                     min_cor = 0,
                                     max_cor = 0.8,
                                     max_iter = 1000,
                                     shrink_step = 0.01) {
  
  if (nConstructs < 2) {
    stop("nConstructs must be >= 2.")
  }
  
  is_pd <- function(M) {
    ev <- eigen(M, symmetric = TRUE, only.values = TRUE)$values
    all(ev > 1e-8)
  }
  
  for (iter in seq_len(max_iter)) {
    
    # Initial random symmetric candidate
    Phi <- diag(1, nConstructs)
    upper_idx <- upper.tri(Phi)
    Phi[upper_idx] <- runif(sum(upper_idx), min = min_cor, max = max_cor)
    Phi[lower.tri(Phi)] <- t(Phi)[lower.tri(Phi)]
    
    # If already valid, return it
    if (is_pd(Phi)) {
      return(Phi)
    }
    
    # Otherwise shrink toward identity until valid
    # alpha = 1   -> original matrix
    # alpha = 0   -> identity matrix
    for (alpha in seq(0.99, 0, by = -shrink_step)) {
      Phi_shrunk <- alpha * Phi + (1 - alpha) * diag(nConstructs)
      
      if (is_pd(Phi_shrunk)) {
        diag(Phi_shrunk) <- 1
        return(Phi_shrunk)
      }
    }
  }
  
  stop("Could not generate a positive definite factor correlation matrix even after shrinkage.")
}

# ------------------------------------------------------------
# Helper: convert continuous data to 7-point Likert
# Uses fixed standard-normal cutpoints after z-standardizing
# each item separately.
# ------------------------------------------------------------
likert7_fixed_mat <- function(dat) {
  
  if (!is.data.frame(dat) && !is.matrix(dat)) {
    stop("dat must be a data.frame or matrix.")
  }
  
  dat <- as.data.frame(dat)
  
  # 7 categories -> 8 cutpoints
  cuts <- qnorm(seq(0, 1, length.out = 8))
  
  out <- lapply(dat, function(x) {
    x <- as.numeric(x)
    
    # Guard against zero variance, just in case
    if (sd(x, na.rm = TRUE) == 0 || is.na(sd(x, na.rm = TRUE))) {
      return(rep(4L, length(x)))
    }
    
    z <- as.numeric(scale(x))
    as.integer(cut(z, breaks = cuts, labels = FALSE, include.lowest = TRUE))
  })
  
  out <- as.data.frame(out)
  names(out) <- names(dat)
  out
}

# ------------------------------------------------------------
# Main function
# Simulates "good respondents" answering a possibly imperfect
# questionnaire:
# - weak/strong loadings allowed
# - factor redundancy allowed via high factor correlations
# - no careless responding / no weird response styles
# ------------------------------------------------------------
simulated_good_responses <- function(nConstructs, nItems, n = 1000, seed = NULL) {
  
  if (!is.null(seed)) set.seed(seed)
  
  # --------------------------
  # Input checks
  # --------------------------
  if (length(nConstructs) != 1 || !is.numeric(nConstructs) || nConstructs < 2 || nConstructs %% 1 != 0) {
    stop("nConstructs must be a single integer >= 2.")
  }
  
  if (!is.numeric(nItems) || length(nItems) != nConstructs || any(nItems < 3) || any(nItems %% 1 != 0)) {
    stop("nItems must be an integer vector of length nConstructs, with all values >= 3.")
  }
  
  if (length(n) != 1 || !is.numeric(n) || n < 1 || n %% 1 != 0) {
    stop("n must be a single integer >= 1.")
  }
  
  # --------------------------
  # Names
  # --------------------------
  factor_names <- paste0("f", seq_len(nConstructs))
  item_names   <- paste0("y", seq_len(sum(nItems)))
  
  # --------------------------
  # Random loadings
  # Range: 0.2–0.8  (shared variance 4%–64%)
  # Rationale: includes weak items that are common in real questionnaires.
  #   The longstring bug fix (2026-03-10e) resolved the issue that previously
  #   made low loadings problematic for detection.
  # --------------------------
  loadings_vec <- round(runif(sum(nItems), min = 0.2, max = 0.8), 2)
  
  # Residual variances so that Var(y) ~= 1 when Var(f)=1
  resid_var_vec <- round(1 - loadings_vec^2, 4)
  
  # --------------------------
  # Random factor correlation matrix
  # Assumption: from 0 to 0.8
  # --------------------------
  Phi <- random_factor_cor_matrix(
    nConstructs = nConstructs,
    min_cor = 0,
    max_cor = 0.8
  )
  
  # --------------------------
  # Build lavaan model syntax
  # --------------------------
  model_lines <- character(0)
  
  current_item <- 1
  loading_table <- data.frame(
    item = character(0),
    factor = character(0),
    loading = numeric(0),
    resid_var = numeric(0),
    stringsAsFactors = FALSE
  )
  
  # Factor measurement part
  for (i in seq_len(nConstructs)) {
    these_items <- item_names[current_item:(current_item + nItems[i] - 1)]
    these_loads <- loadings_vec[current_item:(current_item + nItems[i] - 1)]
    
    line <- paste0(
      factor_names[i],
      " =~ ",
      paste0(these_loads, "*", these_items, collapse = " + ")
    )
    
    model_lines <- c(model_lines, line)
    
    loading_table <- rbind(
      loading_table,
      data.frame(
        item = these_items,
        factor = factor_names[i],
        loading = these_loads,
        resid_var = resid_var_vec[current_item:(current_item + nItems[i] - 1)],
        stringsAsFactors = FALSE
      )
    )
    
    current_item <- current_item + nItems[i]
  }
  
  # Fix latent variances to 1
  model_lines <- c(
    model_lines,
    paste0(factor_names, " ~~ 1*", factor_names)
  )
  
  # Residual variances
  model_lines <- c(
    model_lines,
    paste0(item_names, " ~~ ", resid_var_vec, "*", item_names)
  )
  
  # Factor correlations
  comb_cos <- combn(seq_len(nConstructs), 2)
  for (k in seq_len(ncol(comb_cos))) {
    i <- comb_cos[1, k]
    j <- comb_cos[2, k]
    
    model_lines <- c(
      model_lines,
      paste0(factor_names[i], " ~~ ", round(Phi[i, j], 2), "*", factor_names[j])
    )
  }
  
  model <- paste(model_lines, collapse = "\n")
  
  # --------------------------
  # Simulate continuous latent-response data
  # --------------------------
  dat_cont <- lavaan::simulateData(
    model = model,
    sample.nobs = n,
    meanstructure = FALSE
  )
  
  # --------------------------
  # Convert to 7-point Likert
  # --------------------------
  dat_likert <- likert7_fixed_mat(dat_cont)
  
  # --------------------------
  # Attach useful metadata
  # --------------------------
  attr(dat_likert, "model_syntax") <- model
  attr(dat_likert, "loading_table") <- loading_table
  attr(dat_likert, "factor_cor_matrix") <- Phi
  
  return(dat_likert)
}

