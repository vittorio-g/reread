# Aggiornamenti per Marcello — 2026-04-01

## Decisione architetturale: due versioni del ReReReRe

A partire da oggi il codice espone **due funzioni** nel file `ReReReRe.R`:

| Funzione | Metodo | Default | Status |
|----------|--------|---------|--------|
| `ReReReRe()` | **Standard a 2 livelli**: weighted (tutte le coppie, pesate per \|r\|) per ≤60 item, coupled (top-k% per \|r\|) per >60 item | Sì | **Validato** su dati reali |
| `ReReReRe_F()` | **Per-factor**: EFA → coppie within-factor → pesate per \|r\| | No | **Sperimentale** |

### Perché due versioni?

1. **Standard (`ReReReRe`)** vince sui dataset di validazione esterna con ground truth (AUC=0.635 vs EFA-D=0.561, coupled vince 5/6 dataset)
2. **Per-factor (`ReReReRe_F`)** vince su dati simulati (MCC=0.333 vs 0.305 overall) e su questionari corti reali (es. Pennycook)
3. La discrepanza simulazione→reale dipende dalla qualità dell'EFA su dati reali con cross-loadings e strutture complesse

---

## Metodi EFA testati e scartati come default

### Opzione A — EFA + loading weights (λ_i × λ_j)

54 condizioni simulate (nF=4-20, ipf=3-10, 3 reps). Risultato: vince solo a 60-100 item (+0.027 MCC). I prodotti dei loading sono instabili quando la parallel analysis sottostima nF (~70% del vero).

**Scartata.**

### Opzione D (= EFA-D) — EFA within-factor + |r| weights

360 condizioni simulate (8 nF × 3 ipf × 3 pct × 5 reps). Risultati su dati simulati:

| Range item | Standard | Weighted | EFA-D | Vincitore |
|-----------|----------|----------|-------|-----------|
| <30 | 0.108 | **0.160** | 0.158 | Weighted |
| 30-60 | 0.198 | 0.244 | **0.261** | EFA-D |
| 60-100 | 0.306 | 0.307 | **0.351** | EFA-D |
| 100-200 | 0.501 | 0.405 | **0.529** | EFA-D |
| >200 | **0.656** | 0.393 | 0.653 | Standard |

Overall simulazione: **EFA-D=0.347 > Std=0.303 > Wt=0.286**. Ma su dati reali:

| Dataset | Items | EFA-D AUC | Coupled AUC | EFA-D MCC₁.₅ | Coupled MCC₁.₅ |
|---------|-------|-----------|-------------|-------------|----------------|
| Schroeders | 60 | 0.583 | **0.613** | 0.126 | **0.182** |
| Schneider | 31 | 0.459 | **0.729** | -0.019 | **0.114** |
| Niessen | 100 | 0.587 | **0.637** | 0.037 | **0.073** |
| Goldammer S1 | 60 | 0.605 | **0.717** | 0.167 | **0.314** |
| Goldammer S2 | 60 | 0.599 | **0.656** | 0.096 | **0.286** |
| Goldammer S3 | 60 | **0.532** | 0.456 | 0.096 | -0.032 |

**Media: EFA-D AUC=0.561 vs Coupled AUC=0.635.** Il coupled vince 5/6 dataset.

**Scartata come default.** Rimane disponibile via `ReReReRe_F()` o `mode="efa_d"`.

### Per-factor coherence (Idea 1 da buone_idee.md)

Calcola z-score per ciascun fattore EFA, poi aggrega (mean_z, min_z, combined, prop_low).

Simulazione (360 condizioni): **PF mean(z)=0.333 > Standard=0.305**, vince sotto i 100 item. Ma su dati reali:

| Dataset | Coupled Oracle MCC | PF_mean Oracle MCC |
|---------|-------------------|--------------------|
| Media 6 dataset | **0.238** | 0.124 |

**Stessa dinamica:** vince su simulato, perde su reale. Per lo stesso motivo (EFA scadente su dati reali).

---

## Test su Pennycook & Rand 2019

Applicati i tre metodi (coupled, weighted, EFA-D) ai dati del paper "Lazy, not biased" (Study 1: 30 item, N=782; Study 2: 24 item, N=2564). Flagging capped al 20%.

### Studio 1 — r(CRT, Discernment): paper = .27

| Metodo | Δr | ΔR² | Risultato |
|--------|-----|-----|-----------|
| **EFA-D** | **+0.047** (r .268→.316) | **+0.037** | **Migliore** |
| Weighted | +0.024 (r .268→.293) | +0.024 | Buono |
| Coupled | −0.009 | −0.002 | Peggiora! |

### Studio 2 — r(CRT, Discernment): paper = .21

| Metodo | Δr | ΔR² |
|--------|-----|-----|
| Weighted | +0.019 | +0.012 |
| **EFA-D** | **+0.017** | **+0.010** |
| Coupled | −0.002 | −0.004 |

### Correlazioni Table 1 (Studio 1, 20% flagging)

|  | L_Disc | C_Disc | N_Disc |
|--|--------|--------|--------|
| **Baseline Clinton** | .286 | .090 | .188 |
| **EFA-D Clinton** | **.315** | **.115** | **.231** |
| **Baseline Trump** | .225 | .197 | .232 |
| **EFA-D Trump** | **.301** | **.270** | .244 |

Tutte le correlazioni migliorano. Il guadagno maggiore è per Trump Rep-consistent discernment (.197→.270).

### Conclusione Pennycook

**Su questionari corti, EFA-D e weighted funzionano e migliorano i risultati. Coupled non funziona.** Questo è coerente con le simulazioni: sotto 60 item il coupled non ha abbastanza coppie ad alta |r| per produrre un segnale pulito.

---

## Nota sulla qualità dei dataset di validazione esterna

I dataset con ground truth di careless responding presentano problemi:
- **Goldammer S1-S3**: 63-67% careless (sperimentale, instructed) — rate irrealisticamente alto che corrompe la matrice di correlazione e degrada l'EFA
- **Niessen**: ground truth = speed manipulation, non necessariamente inconsistenza
- **Schneider**: ground truth = classe latente (algoritmica, non sperimentale)
- **Schroeders**: il più pulito, ma 40% careless instructed

La differenza tra "vince su simulazione" e "perde su dati reali" potrebbe essere almeno in parte dovuta alla qualità discutibile del ground truth, non solo ai limiti dell'EFA. Per questo abbiamo aggiunto il test Pennycook: lì il criterio è "i risultati del paper migliorano?", che è un test più ecologicamente valido.

---

## Raccomandazioni pratiche (invariate)

- **corProp = 0.03** (default)
- **z_threshold = 1.5** (robusto universale)
- **auto_z = TRUE** per comodità (la funzione sceglie z basandosi su total_items)
- Usare `ReReReRe()` per analisi primaria
- Usare `ReReReRe_F()` come analisi complementare, specialmente per questionari corti
- Sweet spot: ≥120 item, n ≥ 300

---

## File e struttura

### File core (root)

| File | Funzione | Scopo |
|------|----------|-------|
| `ReReReRe.R` | `ReReReRe()` + `ReReReRe_F()` | Algoritmo principale (standard + per-factor) |
| `Synthetic_Good_Responses_2.R` | `simulated_good_responses()` | Genera dati questionario puliti (CFA-based) |
| `Careless_machine_2.R` | `inject_careless()` | Inietta careless responding |

### Script di test (root)

| File | Scopo |
|------|-------|
| `Test_PerFactor_v2.R` | 360 condizioni per-factor vs standard su dati simulati |
| `Test_PerFactor_RealData.R` | Per-factor vs standard su 6 dataset reali |
| `Test_EFA_D_Full.R` | 360 condizioni EFA-D vs standard/weighted |
| `Test_EFA_D_RealData_NoJohnson.R` | EFA-D vs standard su 6 dataset reali |
| `Test_EFA_OptionD.R` | Confronto 4 metodi (54 condizioni) |

### Dataset

| Cartella | Dataset | Scopo |
|----------|---------|-------|
| `Dataset/dataset_1/` | Pennycook & Rand 2019 | Test applicativo (i risultati migliorano?) |
| `external_datasets/` | 6 dataset + Johnson | Validazione con ground truth |
