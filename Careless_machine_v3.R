#' Inject careless responses — v3 (unified semantics across patterns)
#'
#' Changes from v2 (Simulation_VariancePenalty.R inline version):
#'   - pure_straight and acquiescent now RESPECT corruption_level (were always 100%)
#'   - acquiescent is now SYSTEMATIC: picks one fixed high value, applies +/-1 noise
#'     (was: random sample from high levels, producing sd ~ 0.5)
#'   - Added "fatigue" pattern: corruption concentrated at the end of the item list
#'   - Pattern semantics unified: corruption_level controls how much of the row
#'     is replaced; the pattern controls HOW it is replaced
#'
#' Patterns (all respect corruption_level):
#'   - random         : k random items -> random Likert values
#'   - longstring     : k random items -> 2-5 chunks, each chunk one constant value
#'                      (a.k.a. "scattered-longstring", realistic for scrambled surveys)
#'   - pure_straight  : k random items -> ONE constant value (all same)
#'   - acquiescent    : k random items -> one fixed high value with +/-1 noise
#'                      (systematic "always agree" bias)
#'   - mixed          : half-longstring, half-random on k random items
#'   - fatigue        : LAST k items replaced (position-dependent), mechanical pattern
#'                      resembling longstring (2-5 chunks) at the end

inject_careless_v3 <- function(data,
                                pct_careless = 0.20,
                                pct_types = c(random        = 1/6,
                                              longstring    = 1/6,
                                              pure_straight = 1/6,
                                              acquiescent   = 1/6,
                                              mixed         = 1/6,
                                              fatigue       = 1/6),
                                careless_levels = seq(0.1, 1.0, 0.1),
                                seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  X <- as.data.frame(data); N <- nrow(X); J <- ncol(X)
  for (jj in seq_len(J)) if (!is.numeric(X[[jj]]))
    X[[jj]] <- as.numeric(as.character(X[[jj]]))
  vals_all <- unlist(X); vals_all <- vals_all[!is.na(vals_all)]
  LEVELS <- sort(unique(as.integer(round(vals_all))))
  if (length(LEVELS) < 2 || length(LEVELS) > 15)
    LEVELS <- seq(min(LEVELS), max(LEVELS))
  high_LEVELS <- LEVELS[LEVELS >= quantile(LEVELS, 0.65)]
  high_median <- if (length(high_LEVELS) > 0) median(high_LEVELS) else max(LEVELS)

  pct_types <- pct_types / sum(pct_types)
  M <- round(N * pct_careless)
  type_names <- names(pct_types)
  counts <- setNames(floor(M * pct_types), type_names)
  rem <- M - sum(counts)
  if (rem > 0) counts[seq_len(rem)] <- counts[seq_len(rem)] + 1

  if (M == 0) {
    labels <- data.frame(careless = rep(FALSE, N), careless_pct = 0,
                         pattern = "clean", stringsAsFactors = FALSE)
    return(list(data_corrupted = X, labels = labels))
  }

  careless_ids <- sample(seq_len(N), M, replace = FALSE)

  pat_vec <- character(0); pct_vec <- numeric(0)
  for (pat in type_names) {
    Mp <- counts[pat]
    if (Mp == 0) next
    L <- length(careless_levels)
    cntL <- rep(floor(Mp/L), L); rem2 <- Mp - sum(cntL)
    if (rem2 > 0) cntL[sample(L, rem2)] <- cntL[sample(L, rem2)] + 1
    pcts <- unlist(mapply(function(p, n) rep(p, n), careless_levels, cntL,
                          SIMPLIFY = FALSE))
    if (length(pcts) > 1) pcts <- sample(pcts)
    pat_vec <- c(pat_vec, rep(pat, Mp)); pct_vec <- c(pct_vec, pcts)
  }
  perm <- sample(seq_along(pat_vec))
  pat_vec <- pat_vec[perm]; pct_vec <- pct_vec[perm]

  # ================== Corruption functions ==================
  do_random <- function(row, k) {
    if (k == 0) return(list(row = row, idx = integer(0)))
    idx <- sample(J, k)
    row[idx] <- sample(LEVELS, k, replace = TRUE)
    list(row = row, idx = sort(idx))
  }

  do_longstring <- function(row, k) {
    if (k == 0) return(list(row = row, idx = integer(0)))
    idx <- sample(J, k)
    n_chunks <- if (k >= 4) sample(2:min(5, k), 1) else if (k >= 2) 2 else 1
    chunk <- sample(rep(seq_len(n_chunks), length.out = k))
    for (ch in seq_len(n_chunks)) {
      row[idx[chunk == ch]] <- sample(LEVELS, 1)
    }
    list(row = row, idx = sort(idx))
  }

  # FIXED: pure_straight now respects corruption_level
  # Replace k random items with ONE constant value (same across all k)
  do_pure_straight <- function(row, k) {
    if (k == 0) return(list(row = row, idx = integer(0)))
    idx <- sample(J, k)
    val <- sample(LEVELS, 1)
    row[idx] <- val
    list(row = row, idx = sort(idx))
  }

  # FIXED: acquiescent now respects corruption_level AND is systematic
  # Pick ONE fixed high value (e.g., 6 on 7-point); apply +/-1 integer noise
  do_acquiescent <- function(row, k) {
    if (k == 0) return(list(row = row, idx = integer(0)))
    idx <- sample(J, k)
    base_val <- sample(high_LEVELS, 1)
    noise <- sample(c(-1, 0, 0, 0, 1), k, replace = TRUE)  # mostly no noise
    vals <- pmin(max(LEVELS), pmax(min(LEVELS), base_val + noise))
    row[idx] <- vals
    list(row = row, idx = sort(idx))
  }

  do_mixed <- function(row, k) {
    if (k <= 1) return(do_random(row, k))
    k_long <- max(1, round(k/2)); k_rand <- k - k_long
    L <- do_longstring(row, k_long); row2 <- L$row
    avail <- setdiff(seq_len(J), L$idx)
    if (length(avail) > 0 && k_rand > 0) {
      ri <- sample(avail, min(k_rand, length(avail)))
      row2[ri] <- sample(LEVELS, length(ri), replace = TRUE)
      return(list(row = row2, idx = sort(unique(c(L$idx, ri)))))
    }
    list(row = row2, idx = L$idx)
  }

  # NEW: fatigue = last k items replaced with longstring-like pattern (end-loaded)
  do_fatigue <- function(row, k) {
    if (k == 0) return(list(row = row, idx = integer(0)))
    idx <- (J - k + 1):J  # last k positions
    n_chunks <- if (k >= 4) sample(2:min(5, k), 1) else if (k >= 2) 2 else 1
    chunk <- rep(seq_len(n_chunks), length.out = k)  # sequential chunks
    for (ch in seq_len(n_chunks)) {
      row[idx[chunk == ch]] <- sample(LEVELS, 1)
    }
    list(row = row, idx = sort(idx))
  }

  Xc <- X
  labels <- data.frame(careless = rep(FALSE, N), careless_pct = 0,
                       pattern = "clean", stringsAsFactors = FALSE)
  for (m in seq_len(M)) {
    i <- careless_ids[m]; pct_i <- pct_vec[m]; pat_i <- pat_vec[m]
    k <- min(J, max(0, round(J * pct_i)))
    orig <- as.numeric(Xc[i, , drop = TRUE])
    res <- switch(pat_i,
      "random"        = do_random(orig, k),
      "longstring"    = do_longstring(orig, k),
      "pure_straight" = do_pure_straight(orig, k),
      "acquiescent"   = do_acquiescent(orig, k),
      "mixed"         = do_mixed(orig, k),
      "fatigue"       = do_fatigue(orig, k))
    Xc[i, ] <- res$row
    labels$careless[i] <- TRUE
    labels$careless_pct[i] <- length(res$idx) / J
    labels$pattern[i] <- pat_i
  }
  list(data_corrupted = Xc, labels = labels)
}
