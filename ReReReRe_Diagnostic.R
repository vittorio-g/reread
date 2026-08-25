#' ReReReRe_Diagnostic.R
#'
#' Diagnostico struttura: misura quanto un questionario ha "abbastanza
#' segnale multi-costrutto" perchè ReReReRe possa funzionare bene.
#'
#' Idea: ReReReRe seleziona le top-k% coppie a |r| alto. Se la struttura
#' è ricca, il |r| medio delle top-k% e' molto piu' alto della media |r|
#' globale (rapporto >> 1). Se il questionario non ha struttura
#' multi-costrutto, il rapporto e' ~ 1 e ReReReRe sta solo "rincorrendo
#' il rumore" delle correlazioni casuali.
#'
#' Restituisce:
#'   - separation_ratio : mean(|r|_top_k) / mean(|r|_all)   (>1 buono)
#'   - first_eig_ratio  : eigenvalue_1(top_k_items) / mean(eigenvalues)
#'                        (>2 buono, <1.5 segnale debole)
#'   - n_top_k_items    : numero item unici toccati dalle top-k% coppie
#'   - rec              : raccomandazione testuale
#'   - level            : "ok" / "marginal" / "weak" / "no_structure"

diagnose_structure <- function(data, corProp = 0.03,
                               align_signs = TRUE,
                               warn = TRUE) {
  X <- as.matrix(data)
  J <- ncol(X)
  if (J < 6) {
    return(list(level = "weak", n_pairs = 0, separation_ratio = NA,
                first_eig_ratio = NA, n_top_k_items = J,
                rec = "Questionario troppo corto (<6 item) per il diagnostico."))
  }

  # Matrice di correlazione (con eventuale align_signs trasformato)
  Xn <- X
  if (align_signs) {
    # Senza saperne le keying directions, semplice trick: prendi la median
    # signed correlation di ogni item rispetto al "centroide" e flippa quelli
    # negativi (pseudo-reverse-coding automatico).
    R0 <- suppressWarnings(cor(Xn, use = "pairwise.complete.obs"))
    R0[is.na(R0)] <- 0
    item_med <- apply(R0, 2, function(v) median(v[abs(v) > 0]))
    flip_idx <- which(item_med < 0)
    if (length(flip_idx) > 0) {
      m_max <- apply(Xn[, flip_idx, drop = FALSE], 2,
                     function(v) max(v, na.rm = TRUE))
      m_min <- apply(Xn[, flip_idx, drop = FALSE], 2,
                     function(v) min(v, na.rm = TRUE))
      for (k in seq_along(flip_idx)) {
        Xn[, flip_idx[k]] <- (m_max[k] + m_min[k]) - Xn[, flip_idx[k]]
      }
    }
  }

  R <- suppressWarnings(cor(Xn, use = "pairwise.complete.obs"))
  R[is.na(R)] <- 0
  diag(R) <- NA
  abs_r <- abs(R[upper.tri(R)])
  abs_r <- abs_r[!is.na(abs_r)]
  if (length(abs_r) == 0) {
    return(list(level = "no_structure", n_pairs = 0,
                separation_ratio = NA, first_eig_ratio = NA,
                n_top_k_items = 0,
                rec = "Nessuna coppia di item con varianza informativa."))
  }

  k_top <- max(1, round(corProp * length(abs_r)))
  threshold <- sort(abs_r, decreasing = TRUE)[k_top]
  top_mean <- mean(abs_r[abs_r >= threshold])
  all_mean <- mean(abs_r)
  separation_ratio <- if (all_mean > 0) top_mean / all_mean else NA

  # First eigenvalue su sottomatrice degli item nelle top-k% coppie
  abs_full <- abs(R)
  abs_full[is.na(abs_full)] <- 0
  diag(abs_full) <- 0
  pair_idx <- which(abs_full >= threshold, arr.ind = TRUE)
  unique_items <- unique(c(pair_idx[, 1], pair_idx[, 2]))
  n_top_k_items <- length(unique_items)

  first_eig_ratio <- NA
  if (n_top_k_items >= 4) {
    R_sub <- R[unique_items, unique_items]
    R_sub[is.na(R_sub)] <- 0
    diag(R_sub) <- 1
    eg <- tryCatch(eigen(R_sub, symmetric = TRUE, only.values = TRUE)$values,
                   error = function(e) NULL)
    if (!is.null(eg)) {
      first_eig_ratio <- eg[1] / mean(eg)
    }
  }

  # Classifica in 4 livelli
  level <- if (is.na(separation_ratio)) "no_structure"
           else if (separation_ratio < 1.5) "no_structure"
           else if (separation_ratio < 2.0) "weak"
           else if (separation_ratio < 3.0) "marginal"
           else "ok"

  rec <- switch(level,
    "ok" = paste0("Struttura adeguata (separation ratio = ",
                  round(separation_ratio, 2), "). ReReReRe funzionera'."),
    "marginal" = paste0("Struttura marginale (separation ratio = ",
                        round(separation_ratio, 2),
                        "). ReReReRe puo' essere usato ma il MCC sara' ridotto. ",
                        "Considera anche IRV+LongString."),
    "weak" = paste0("Struttura DEBOLE (separation ratio = ",
                    round(separation_ratio, 2),
                    "). ReReReRe potrebbe produrre molti falsi positivi. ",
                    "Raccomandato: IRV + LongString come metodo principale, ",
                    "ReReReRe come check secondario."),
    "no_structure" = paste0("Struttura ASSENTE (separation ratio = ",
                            round(separation_ratio, 2),
                            "). ReReReRe NON e' adatto a questo questionario. ",
                            "Usa esclusivamente IRV + LongString.")
  )

  if (warn && level %in% c("weak", "no_structure")) {
    warning(rec, call. = FALSE)
  }

  list(level = level,
       n_items = J,
       n_pairs = length(abs_r),
       n_top_k_pairs = k_top,
       n_top_k_items = n_top_k_items,
       separation_ratio = separation_ratio,
       first_eig_ratio = first_eig_ratio,
       top_k_mean_r = top_mean,
       all_mean_r = all_mean,
       rec = rec)
}
