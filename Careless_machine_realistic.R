#' Careless_machine_realistic.R
#'
#' Estensione "realistic" che aggiunge 3 nuovi pattern di careless al set
#' v3 esistente. I 6 pattern legacy (random, longstring, pure_straight,
#' acquiescent, mixed, fatigue) sono mantenuti identici al
#' Careless_machine_v3.R; ne aggiungo 3 ispirati alla letteratura empirica:
#'
#'   - anchor_noise        : ancora a 3/4/5 + N(0, 0.7), arrotondato a Likert.
#'                           Pattern rumoroso non uniforme (IRV non lo vede).
#'   - language_barrier    : per ogni item, prob. (k/J) di "non capirlo".
#'                           Numero di item corrotti ~ Binomial(J, k/J).
#'   - positional_fatigue  : prob crescente linearmente nella posizione
#'                           (0 all'inizio, 2*k/J alla fine). Modella fatica.
#'
#' Per discriminare dal "fatigue" legacy (che mette gli ultimi k item costanti),
#' positional_fatigue è probabilistico e RUMOROSO (random sostitutivo).

inject_careless_realistic <- function(
    data,
    pct_careless = 0.20,
    pct_types = c(random             = 1/9,
                  longstring         = 1/9,
                  pure_straight      = 1/9,
                  acquiescent        = 1/9,
                  mixed              = 1/9,
                  fatigue            = 1/9,
                  anchor_noise       = 1/9,
                  language_barrier   = 1/9,
                  positional_fatigue = 1/9),
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
  L_min <- min(LEVELS); L_max <- max(LEVELS)
  high_LEVELS <- LEVELS[LEVELS >= quantile(LEVELS, 0.65)]
  high_median <- if (length(high_LEVELS) > 0) median(high_LEVELS) else max(LEVELS)
  midpoint <- median(LEVELS)

  pct_types <- pct_types / sum(pct_types)
  M <- round(N * pct_careless)
  type_names <- names(pct_types)
  counts <- setNames(floor(M * pct_types), type_names)
  rem <- M - sum(counts)
  if (rem > 0) counts[seq_len(rem)] <- counts[seq_len(rem)] + 1

  if (M == 0) {
    return(list(data_corrupted = X,
                labels = data.frame(careless = rep(FALSE, N),
                                    careless_pct = 0,
                                    pattern = "clean",
                                    stringsAsFactors = FALSE)))
  }

  careless_ids <- sample(seq_len(N), M)
  pat_vec <- character(0); pct_vec <- numeric(0)
  for (pat in type_names) {
    Mp <- counts[pat]; if (Mp == 0) next
    L <- length(careless_levels)
    cntL <- rep(floor(Mp / L), L); rem2 <- Mp - sum(cntL)
    if (rem2 > 0) cntL[sample(L, rem2)] <- cntL[sample(L, rem2)] + 1
    pcts <- unlist(mapply(function(p, n) rep(p, n), careless_levels, cntL,
                          SIMPLIFY = FALSE))
    if (length(pcts) > 1) pcts <- sample(pcts)
    pat_vec <- c(pat_vec, rep(pat, Mp)); pct_vec <- c(pct_vec, pcts)
  }
  perm <- sample(seq_along(pat_vec))
  pat_vec <- pat_vec[perm]; pct_vec <- pct_vec[perm]

  # ============ Pattern legacy ============
  do_random <- function(row, k) {
    if (k == 0) return(list(row = row, idx = integer(0)))
    idx <- sample(J, k); row[idx] <- sample(LEVELS, k, replace = TRUE)
    list(row = row, idx = sort(idx))
  }
  do_longstring <- function(row, k) {
    if (k == 0) return(list(row = row, idx = integer(0)))
    idx <- sample(J, k)
    n_chunks <- if (k >= 4) sample(2:min(5, k), 1) else if (k >= 2) 2 else 1
    chunk <- sample(rep(seq_len(n_chunks), length.out = k))
    for (ch in seq_len(n_chunks)) row[idx[chunk == ch]] <- sample(LEVELS, 1)
    list(row = row, idx = sort(idx))
  }
  do_pure_straight <- function(row, k) {
    if (k == 0) return(list(row = row, idx = integer(0)))
    idx <- sample(J, k); row[idx] <- sample(LEVELS, 1)
    list(row = row, idx = sort(idx))
  }
  do_acquiescent <- function(row, k) {
    if (k == 0) return(list(row = row, idx = integer(0)))
    idx <- sample(J, k)
    base_val <- sample(high_LEVELS, 1)
    noise <- sample(c(-1, 0, 0, 0, 1), k, replace = TRUE)
    row[idx] <- pmin(L_max, pmax(L_min, base_val + noise))
    list(row = row, idx = sort(idx))
  }
  do_mixed <- function(row, k) {
    if (k <= 1) return(do_random(row, k))
    k_long <- max(1, round(k / 2)); k_rand <- k - k_long
    L <- do_longstring(row, k_long); row2 <- L$row
    avail <- setdiff(seq_len(J), L$idx)
    if (length(avail) > 0 && k_rand > 0) {
      ri <- sample(avail, min(k_rand, length(avail)))
      row2[ri] <- sample(LEVELS, length(ri), replace = TRUE)
      return(list(row = row2, idx = sort(unique(c(L$idx, ri)))))
    }
    list(row = row2, idx = L$idx)
  }
  do_fatigue <- function(row, k) {
    if (k == 0) return(list(row = row, idx = integer(0)))
    idx <- (J - k + 1):J
    n_chunks <- if (k >= 4) sample(2:min(5, k), 1) else if (k >= 2) 2 else 1
    chunk <- rep(seq_len(n_chunks), length.out = k)
    for (ch in seq_len(n_chunks)) row[idx[chunk == ch]] <- sample(LEVELS, 1)
    list(row = row, idx = sort(idx))
  }

  # ============ Pattern nuovi (realistic) ============
  do_anchor_noise <- function(row, k) {
    if (k == 0) return(list(row = row, idx = integer(0)))
    anchor <- sample(c(midpoint - 1, midpoint, midpoint + 1), 1)
    anchor <- max(L_min, min(L_max, anchor))
    idx <- sample(J, k)
    noisy <- anchor + rnorm(k, 0, 0.7)
    row[idx] <- as.integer(pmax(L_min, pmin(L_max, round(noisy))))
    list(row = row, idx = sort(idx))
  }
  do_language_barrier <- function(row, k) {
    if (k == 0) return(list(row = row, idx = integer(0)))
    p <- k / J
    idx <- which(runif(J) < p)
    if (length(idx) == 0) idx <- sample(J, 1)
    row[idx] <- sample(LEVELS, length(idx), replace = TRUE)
    list(row = row, idx = sort(idx))
  }
  do_positional_fatigue <- function(row, k) {
    if (k == 0) return(list(row = row, idx = integer(0)))
    base_p <- k / J
    pos_p  <- pmin(1, 2 * base_p * (seq_len(J) / J))
    idx <- which(runif(J) < pos_p)
    if (length(idx) == 0) idx <- J
    row[idx] <- sample(LEVELS, length(idx), replace = TRUE)
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
      "random"             = do_random(orig, k),
      "longstring"         = do_longstring(orig, k),
      "pure_straight"      = do_pure_straight(orig, k),
      "acquiescent"        = do_acquiescent(orig, k),
      "mixed"              = do_mixed(orig, k),
      "fatigue"            = do_fatigue(orig, k),
      "anchor_noise"       = do_anchor_noise(orig, k),
      "language_barrier"   = do_language_barrier(orig, k),
      "positional_fatigue" = do_positional_fatigue(orig, k))
    Xc[i, ] <- res$row
    labels$careless[i] <- TRUE
    labels$careless_pct[i] <- length(res$idx) / J
    labels$pattern[i] <- pat_i
  }
  list(data_corrupted = Xc, labels = labels)
}
