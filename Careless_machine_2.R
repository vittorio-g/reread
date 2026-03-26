#' Inietta risposte careless in un dataset item-level (versione con livelli custom)
#'
#' @param data data.frame/matrice (righe = rispondenti, colonne = item)
#' @param pct_careless percentuale di rispondenti careless (0-1 o 0-100)
#' @param pct_types percentuali per c("random","longstring","mixed") (0-1 o 0-100)
#' @param careless_levels vettore dei livelli di corruzione item (% degli item corrotti
#'        per ciascun gruppo, 0-1 o 0-100). Es: c(40,60,90) o c(.4,.6,.9).
#'        Distribuiti equamente *per pattern*.
#' @param seed intero opzionale per riproducibilità
#' @return list(data_corrupted, labels, meta)
inject_careless <- function(data,
                            pct_careless,
                            pct_types = c(random = 1/3, longstring = 1/3, mixed = 1/3),
                            careless_levels = seq(0.1,1,0.1),
                            seed = NULL) {
  
  if (!is.null(seed)) set.seed(seed)
  
  # --- Preparazione e controlli ---
  X <- as.data.frame(data, stringsAsFactors = FALSE)
  N <- nrow(X); J <- ncol(X)
  if (N < 1 || J < 2) stop("Dataset troppo piccolo: servono >=1 riga e >=2 colonne.")
  
  # (Robustezza) prova a convertire colonne non-numeriche in numeriche se sono codificate come stringhe/fattori
  for (jj in seq_len(J)) {
    if (!is.numeric(X[[jj]])) {
      tmp <- suppressWarnings(as.numeric(as.character(X[[jj]])))
      non_na <- !is.na(X[[jj]])
      if (any(non_na)) {
        ok_rate <- sum(!is.na(tmp[non_na])) / sum(non_na)
        if (ok_rate < 0.95) {
          stop(sprintf(
            "Colonna %s non è numerica (o non convertibile in modo affidabile).",
            if (!is.null(names(X)[jj]) && nzchar(names(X)[jj])) names(X)[jj] else jj
          ))
        }
      }
      X[[jj]] <- tmp
    }
  }
  
  # Percentuale careless (0-1 o 0-100)
  p_car <- if (pct_careless > 1) pct_careless/100 else pct_careless
  if (p_car < 0 || p_car > 1) stop("pct_careless deve essere in [0,1] o [0,100].")
  
  # Percentuali per tipo (0-1 o 0-100, poi normalizzate) - FIX: preserva/usa i nomi
  type_names <- c("random","longstring","mixed")
  pct_types <- unlist(pct_types, use.names = TRUE)
  
  if (is.null(names(pct_types))) {
    if (length(pct_types) != 3) {
      stop("pct_types deve avere 3 valori: c(random=..., longstring=..., mixed=...).")
    }
    pct_types <- setNames(as.numeric(pct_types), type_names)
  } else {
    # riordina e verifica che ci siano i 3 nomi richiesti
    if (!all(type_names %in% names(pct_types))) {
      stop("pct_types deve contenere i nomi: random, longstring, mixed.")
    }
    pct_types <- as.numeric(pct_types[type_names])
    names(pct_types) <- type_names
  }
  
  if (any(is.na(pct_types))) stop("pct_types contiene NA non validi.")
  if (any(pct_types < 0)) stop("pct_types non può contenere valori negativi.")
  if (max(pct_types) > 1) pct_types <- pct_types/100
  if (sum(pct_types) <= 0) stop("Somma di pct_types non valida.")
  props <- pct_types / sum(pct_types)  # già nell'ordine type_names
  
  # Livelli di corruzione (0-1 o 0-100) – mantieni l'ordine fornito dall'utente
  levels_pct <- if (max(careless_levels) > 1) careless_levels/100 else careless_levels
  if (any(levels_pct <= 0) || any(levels_pct > 1)) {
    stop("careless_levels deve essere in (0,1] o (0,100]. Esempio: c(40,60,90).")
  }
  
  # Numero rispondenti careless
  M <- round(N * p_car)
  M <- min(M, N)
  if (M == 0) {
    labels <- data.frame(
      careless = rep(FALSE, N),
      careless_pct = rep(0, N),
      pattern = rep("clean", N),
      stringsAsFactors = FALSE
    )
    return(list(
      data_corrupted = X,
      labels = labels,
      meta = list(
        N = N, J = J, M = M,
        counts = c(random=0,longstring=0,mixed=0),
        levels = levels_pct,
        seed = seed,
        levels_used = NA
      )
    ))
  }
  
  # Distribuzione per pattern tra i M careless (gestione arrotondamenti)
  base_counts <- floor(M * props)
  remainder <- M - sum(base_counts)
  if (remainder > 0) {
    order_idx <- order(props, decreasing = TRUE)
    for (k in seq_len(remainder)) {
      ii <- order_idx[((k - 1) %% 3) + 1]
      base_counts[ii] <- base_counts[ii] + 1
    }
  }
  names(base_counts) <- type_names
  
  # Helper: distribuisce Mp partecipanti in L livelli in modo (quasi) equo
  distribute_levels <- function(Mp, L) {
    if (Mp == 0) return(rep(0, L))
    base <- rep(floor(Mp/L), L)
    rem  <- Mp - sum(base)
    if (rem > 0) {
      idx <- sample(seq_len(L), rem)
      base[idx] <- base[idx] + 1
    }
    base
  }
  
  careless_ids <- sample(seq_len(N), M, replace = FALSE)
  
  pattern_vec <- character(0)
  assigned_pcts <- numeric(0)
  per_pattern_level_counts <- list()
  
  for (pat in type_names) {
    Mp <- base_counts[pat]
    if (Mp > 0) {
      L <- length(levels_pct)
      cntL <- distribute_levels(Mp, L)
      per_pattern_level_counts[[pat]] <- setNames(cntL, paste0(round(levels_pct*100, 3), "%"))
      
      pcts_pat <- unlist(mapply(function(pct, nrep) rep(pct, nrep),
                                pct = levels_pct, nrep = cntL, SIMPLIFY = FALSE))
      if (length(pcts_pat) > 1) pcts_pat <- sample(pcts_pat, length(pcts_pat))
      
      pattern_vec   <- c(pattern_vec, rep(pat, Mp))
      assigned_pcts <- c(assigned_pcts, pcts_pat)
    } else {
      per_pattern_level_counts[[pat]] <- setNames(rep(0, length(levels_pct)),
                                                  paste0(round(levels_pct*100, 3), "%"))
    }
  }
  
  if (length(pattern_vec) > 1) {
    perm <- sample(seq_along(pattern_vec))
    pattern_vec   <- pattern_vec[perm]
    assigned_pcts <- assigned_pcts[perm]
  }
  
  # --- Utility: livelli possibili di risposta (Likert globale) ---
  vals_all <- as.integer(round(unlist(X)))
  vals_all <- vals_all[!is.na(vals_all)]
  
  if (length(vals_all) == 0) {
    stop("Impossibile inferire i livelli di risposta: tutte le celle sono NA.")
  }
  
  uniq_vals <- sort(unique(vals_all))
  if (length(uniq_vals) >= 2 && length(uniq_vals) <= 15) {
    LEVELS <- uniq_vals
  } else {
    rng <- range(vals_all, na.rm = TRUE)
    LEVELS <- seq(from = floor(rng[1]), to = ceiling(rng[2]))
  }
  
  # --- Helper: RANDOM ---
  corrupt_random <- function(row_vec, k) {
    if (k <= 0) return(list(row=row_vec, idx=integer(0)))
    idx <- sample(seq_along(row_vec), k, replace = FALSE)
    row_vec[idx] <- sample(LEVELS, k, replace = TRUE)
    list(row=row_vec, idx=idx)
  }
  
  # --- Helper: LONGSTRING (scattered positions, same-value chunks) ---
  # Simulates straight-line responding in a questionnaire with scrambled item order.
  # Selects k RANDOM (non-contiguous) item positions and assigns them to 2-5 chunks,

  # each chunk getting a single constant Likert value. This avoids the confound where
  # contiguous blocks align with factor-ordered items in simulated data.
  corrupt_longstring <- function(row_vec, k) {
    if (k <= 0) return(list(row=row_vec, idx=integer(0)))

    # Select k random (non-contiguous) positions
    idx_all <- sample(seq_along(row_vec), k, replace = FALSE)

    # Decide number of same-value chunks (each chunk = one constant response)
    if (k >= 4) {
      n_chunks <- sample(2:min(5, k), 1)
    } else if (k >= 2) {
      n_chunks <- 2
    } else {
      n_chunks <- 1
    }

    # Assign positions to chunks randomly
    chunk_membership <- sample(rep(seq_len(n_chunks), length.out = k))

    for (ch in seq_len(n_chunks)) {
      idx_chunk <- idx_all[chunk_membership == ch]
      val <- sample(LEVELS, 1)
      row_vec[idx_chunk] <- val
    }

    list(row = row_vec, idx = sort(idx_all))
  }
  
  # --- Helper: MIXED (50% longstring, 50% random) ---
  corrupt_mixed <- function(row_vec, k) {
    if (k <= 1) return(corrupt_random(row_vec, k))
    k_long <- max(1, round(0.5 * k))
    k_rand <- k - k_long
    Ls <- corrupt_longstring(row_vec, k_long)
    row_vec2 <- Ls$row
    used_idx <- Ls$idx
    
    available <- setdiff(seq_along(row_vec2), used_idx)
    if (length(available) > 0 && k_rand > 0) {
      idx_rand <- sample(available, min(k_rand, length(available)), replace = FALSE)
      row_vec2[idx_rand] <- sample(LEVELS, length(idx_rand), replace = TRUE)
      used_idx <- sort(unique(c(used_idx, idx_rand)))
    }
    list(row=row_vec2, idx=used_idx)
  }
  
  # --- Applica la corruzione ---
  Xc <- X
  labels <- data.frame(
    careless = rep(FALSE, N),
    careless_pct = rep(0, N),
    pattern = rep("clean", N),
    stringsAsFactors = FALSE
  )
  
  for (m_idx in seq_len(M)) {
    i <- careless_ids[m_idx]
    pct_i <- assigned_pcts[m_idx]
    pat_i <- pattern_vec[m_idx]
    
    k <- round(J * pct_i)
    k <- max(0, min(k, J))
    
    original_row <- as.numeric(Xc[i, , drop = TRUE])
    
    res <- switch(pat_i,
                  "random"     = corrupt_random(original_row, k),
                  "longstring" = corrupt_longstring(original_row, k),
                  "mixed"      = corrupt_mixed(original_row, k))
    
    Xc[i, ] <- res$row
    labels$careless[i] <- TRUE
    labels$careless_pct[i] <- length(res$idx) / J
    labels$pattern[i] <- pat_i
  }
  
  meta <- list(
    N = N, J = J, M = M,
    counts_by_pattern = base_counts,
    counts_by_pattern_levels = per_pattern_level_counts,
    seed = seed,
    levels_used = LEVELS,
    careless_levels = levels_pct
  )
  
  list(data_corrupted = Xc, labels = labels, meta = meta)
}
