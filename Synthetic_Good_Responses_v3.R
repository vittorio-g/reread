#' Synthetic_Good_Responses_v3.R
#'
#' Estensione "realistica" del generatore CFA per testare ReReReRe in condizioni
#' meno idealizzate. Aggiunge:
#'   - cross-loadings su fattori secondari (proporzione e range parametrizzati)
#'   - method-effect generale (un fattore latente che carica positivamente su tutti
#'     gli item con un loading random in [0, method_strength])
#'   - popolazione mista di stili di risposta NON-careless (normal / extreme /
#'     central / acquiescent_attentive)
#'   - skew opzionale per una frazione di item (cutpoint Likert shiftati)
#'
#' L'output è (deliberatamente) un data.frame con la stessa interfaccia di
#' simulated_good_responses() del simulatore canonico.

library(lavaan)
library(dplyr)

# ------------------------------------------------------------
# Helper riusato: fattore di correlazione random PD
# ------------------------------------------------------------
.random_factor_cor_matrix_v3 <- function(nConstructs,
                                         min_cor = 0,
                                         max_cor = 0.8,
                                         max_iter = 1000,
                                         shrink_step = 0.01) {
  if (nConstructs < 2) stop("nConstructs must be >= 2.")
  is_pd <- function(M) all(eigen(M, symmetric = TRUE, only.values = TRUE)$values > 1e-8)
  for (iter in seq_len(max_iter)) {
    Phi <- diag(1, nConstructs)
    upper_idx <- upper.tri(Phi)
    Phi[upper_idx] <- runif(sum(upper_idx), min = min_cor, max = max_cor)
    Phi[lower.tri(Phi)] <- t(Phi)[lower.tri(Phi)]
    if (is_pd(Phi)) return(Phi)
    for (alpha in seq(0.99, 0, by = -shrink_step)) {
      Phi_shrunk <- alpha * Phi + (1 - alpha) * diag(nConstructs)
      if (is_pd(Phi_shrunk)) { diag(Phi_shrunk) <- 1; return(Phi_shrunk) }
    }
  }
  stop("Could not generate a positive definite factor correlation matrix.")
}

# ------------------------------------------------------------
# Discretizzazione Likert con possibile skew per item
# Per default cutpoint simmetrici qnorm(seq(0,1,8)).
# Con skew_prob > 0, una frazione skew_prob di item viene generata con
# cutpoint shiftati: floor effect (shift positivo) o ceiling effect.
# ------------------------------------------------------------
.likert7_v3 <- function(dat, skew_prob = 0, seed_off = 0) {
  cuts_sym  <- qnorm(seq(0, 1, length.out = 8))
  cuts_floor <- qnorm(seq(0,    0.85, length.out = 8))   # 70% rispondenti in cat 1-2
  cuts_ceil  <- qnorm(seq(0.15, 1,    length.out = 8))   # specchio
  out <- as.data.frame(dat)
  for (jj in seq_len(ncol(out))) {
    x <- as.numeric(out[[jj]])
    if (sd(x, na.rm = TRUE) == 0 || is.na(sd(x, na.rm = TRUE))) {
      out[[jj]] <- rep(4L, length(x)); next
    }
    z <- as.numeric(scale(x))
    cuts_use <- cuts_sym
    if (skew_prob > 0 && runif(1) < skew_prob) {
      cuts_use <- if (runif(1) < 0.5) cuts_floor else cuts_ceil
    }
    out[[jj]] <- as.integer(cut(z, breaks = cuts_use, labels = FALSE,
                                include.lowest = TRUE))
    # Pulisci NA (al limite estremo dei cut shiftati)
    bad <- which(is.na(out[[jj]]))
    if (length(bad)) out[[jj]][bad] <- ifelse(z[bad] < cuts_use[1], 1L, 7L)
  }
  out
}

# ------------------------------------------------------------
# Applica stile di risposta NON-careless prima della discretizzazione
# Modifica i z continui:
#   - normal:                identità
#   - extreme:               z * 1.6  (boost varianza intra-respondente)
#   - central:               z * 0.55 (compressione)
#   - acquiescent_attentive: z + 0.7  (shift positivo costante)
# Il vettore di stili viene generato per N rispondenti dato style_mix.
# ------------------------------------------------------------
.apply_style <- function(dat_cont, style_mix) {
  N <- nrow(dat_cont)
  styles <- c("normal", "extreme", "central", "acquiescent_attentive")
  probs <- c(style_mix$normal,
             style_mix$extreme,
             style_mix$central,
             style_mix$acquiescent_attentive)
  if (sum(probs) == 0) return(list(data = dat_cont, style = rep("normal", N)))
  probs <- probs / sum(probs)
  style_vec <- sample(styles, size = N, replace = TRUE, prob = probs)

  out <- as.matrix(dat_cont)
  for (i in seq_len(N)) {
    s <- style_vec[i]
    if (s == "extreme")                 out[i, ] <- out[i, ] * 1.6
    else if (s == "central")             out[i, ] <- out[i, ] * 0.55
    else if (s == "acquiescent_attentive") out[i, ] <- out[i, ] + 0.7
  }
  list(data = as.data.frame(out), style = style_vec)
}

# ------------------------------------------------------------
# Main function — drop-in replacement realistico
# ------------------------------------------------------------
simulated_good_responses_v3 <- function(
    nConstructs, nItems, n = 1000,
    cross_loading_prob = 0,            # frazione di item con un secondario
    cross_loading_range = c(0.15, 0.35),
    method_strength = 0,               # max loading random sul method factor
    skew_prob = 0,                     # frazione di item con cutpoint Likert skew
    style_mix = list(normal = 1,
                     extreme = 0,
                     central = 0,
                     acquiescent_attentive = 0),
    seed = NULL) {

  if (!is.null(seed)) set.seed(seed)

  # --- Validazioni base
  if (length(nConstructs) != 1 || nConstructs < 2)
    stop("nConstructs must be a single integer >= 2.")
  if (length(nItems) != nConstructs || any(nItems < 3))
    stop("nItems must have length nConstructs and all >= 3.")

  factor_names <- paste0("f", seq_len(nConstructs))
  total_items  <- sum(nItems)
  item_names   <- paste0("y", seq_len(total_items))

  # --- Loadings primari random U(0.2, 0.8)
  primary_loads <- round(runif(total_items, 0.2, 0.8), 2)

  # --- Phi tra fattori sostantivi
  Phi <- .random_factor_cor_matrix_v3(nConstructs, 0, 0.8)

  # --- Costruzione factor->items mapping
  item_factor <- integer(total_items)
  cur <- 1
  for (i in seq_len(nConstructs)) {
    item_factor[cur:(cur + nItems[i] - 1)] <- i
    cur <- cur + nItems[i]
  }

  # --- Cross-loadings: per ogni item, con probabilità cross_loading_prob
  #     aggiungi un loading su un fattore *diverso* dal primario
  cross_secondary <- integer(total_items)
  cross_load      <- numeric(total_items)
  if (cross_loading_prob > 0 && nConstructs >= 2) {
    for (j in seq_len(total_items)) {
      if (runif(1) < cross_loading_prob) {
        candidates <- setdiff(seq_len(nConstructs), item_factor[j])
        cross_secondary[j] <- sample(candidates, 1)
        cross_load[j]      <- round(runif(1, cross_loading_range[1],
                                          cross_loading_range[2]), 2)
      }
    }
  }

  # --- Method effect: se method_strength > 0, fattore "m1" con loading random per item
  use_method <- method_strength > 0
  method_loads <- if (use_method)
    round(runif(total_items, 0, method_strength), 2) else
    rep(0, total_items)

  # --- Calcola residui in modo che Var(y) ~= 1
  # Var(y_j) = primary^2 + cross^2 + method^2 + 2*primary*cross*Phi_pq + resid
  # Approssimo trascurando il termine di covarianza tra primario e secondario
  # (è piccolo e Phi è random): resid = max(eps, 1 - sum di quadrati).
  resid_var_vec <- pmax(0.05,
                        round(1 - primary_loads^2 - cross_load^2 - method_loads^2, 4))

  # --- Costruisce sintassi lavaan
  model_lines <- character(0)

  # Misura: per ogni fattore primario raccogli i suoi item
  for (f in seq_len(nConstructs)) {
    items_f <- which(item_factor == f)
    if (length(items_f) == 0) next
    pieces <- paste0(primary_loads[items_f], "*", item_names[items_f])
    model_lines <- c(model_lines,
                     paste0(factor_names[f], " =~ ",
                            paste(pieces, collapse = " + ")))
  }

  # Misura: cross-loadings
  for (f in seq_len(nConstructs)) {
    items_x <- which(cross_secondary == f & cross_load > 0)
    if (length(items_x) == 0) next
    pieces <- paste0(cross_load[items_x], "*", item_names[items_x])
    model_lines <- c(model_lines,
                     paste0(factor_names[f], " =~ ",
                            paste(pieces, collapse = " + ")))
  }

  # Method factor (carica su tutti gli item se attivo)
  if (use_method) {
    pieces <- paste0(method_loads, "*", item_names)
    # Tieni solo loadings >0 per leggibilità; lavaan accetta anche con 0
    keep <- method_loads > 0
    if (any(keep)) {
      model_lines <- c(model_lines,
                       paste0("m1 =~ ",
                              paste(pieces[keep], collapse = " + ")))
      # Method factor scorrelato dai sostantivi (assunto)
      model_lines <- c(model_lines, "m1 ~~ 1*m1")
      for (f in seq_len(nConstructs))
        model_lines <- c(model_lines,
                         paste0("m1 ~~ 0*", factor_names[f]))
    } else {
      use_method <- FALSE  # nessun loading effettivo
    }
  }

  # Varianza dei fattori sostantivi = 1
  model_lines <- c(model_lines, paste0(factor_names, " ~~ 1*", factor_names))

  # Residui item
  model_lines <- c(model_lines,
                   paste0(item_names, " ~~ ", resid_var_vec, "*", item_names))

  # Correlazioni tra fattori sostantivi (dal Phi)
  for (i in seq_len(nConstructs - 1)) {
    for (j in (i + 1):nConstructs) {
      model_lines <- c(model_lines,
                       paste0(factor_names[i], " ~~ ", round(Phi[i, j], 2),
                              "*", factor_names[j]))
    }
  }

  model <- paste(model_lines, collapse = "\n")

  # --- Simula dati continui
  dat_cont <- lavaan::simulateData(model = model, sample.nobs = n,
                                   meanstructure = FALSE)

  # --- Stile di risposta non-careless
  styled <- .apply_style(dat_cont, style_mix)
  dat_cont <- styled$data

  # --- Discretizzazione Likert (con eventuale skew per item)
  dat_likert <- .likert7_v3(dat_cont, skew_prob = skew_prob)

  # --- Metadata
  attr(dat_likert, "model_syntax") <- model
  attr(dat_likert, "factor_cor_matrix") <- Phi
  attr(dat_likert, "primary_loadings") <- primary_loads
  attr(dat_likert, "cross_secondary") <- cross_secondary
  attr(dat_likert, "cross_loadings") <- cross_load
  attr(dat_likert, "method_loadings") <- method_loads
  attr(dat_likert, "response_style") <- styled$style
  attr(dat_likert, "config") <- list(
    cross_loading_prob = cross_loading_prob,
    cross_loading_range = cross_loading_range,
    method_strength = method_strength,
    skew_prob = skew_prob,
    style_mix = style_mix
  )

  return(dat_likert)
}
