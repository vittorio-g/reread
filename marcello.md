# Aggiornamenti per Marcello — 2026-04-21

## Attacco alla corruzione parziale: variance_penalty + split-half

**Problema identificato (2026-04-21):** sotto GT=corruzione>50% + auto-z, i pattern
random/longstring/mixed a 60-80% di corruzione vengono detectati solo al 66-77%.
Il motivo è la **dilution**: su 60% di item corrotti restano ~16% di coppie tutte
pulite che alzano la coherence individuale, portando il z sopra soglia.

**Due nuove mosse implementate:**

### (1) variance_penalty — integrato in ReReReRe()

Formula: `z_adjusted = z_raw - α * exp(-sd_respondent / β)` con α=3, β=0.5 default.

Chi ha varianza interna bassa (straight-liner: sd=0) prende una penalty piena
(−3.0 sullo z). Chi risponde normalmente (sd≥1.5) non è toccato.

**Risultato sui dati simulati (N=12000):** MCC sale da **0.567 → 0.582** (+0.015).
Guadagno concentrato su:
- pure_straight: 0.84 → 0.93 (+0.09)
- acquiescent: 0.84 → 0.93 (+0.09)

Attivabile con `ReReReRe(data, variance_penalty = TRUE)`.

### (2) score_split_half — nuovo file ReReReRe_SplitHalf.R

Splitta le k coppie coupled in 2 metà, calcola z separato per ciascuna metà
(con propria baseline permutativa), ripete B=5 volte con split diversi, e
aggrega.

**Due comportamenti molto diversi:**

| Aggregazione | Sens | Spec | MCC |
|---|:---:|:---:|:---:|
| `"mean"` | 0.89 | 0.73 | **0.589** (best binary) |
| `"min"` | 0.99 | 0.46 | 0.464 |

**Per-pattern sensitivity sui careless "difficili" (corruption >50%):**

| Pattern | std | std+VP | sh_mean | **sh_min** |
|---------|:---:|:---:|:---:|:---:|
| random | 0.66 | 0.68 | 0.79 | **0.99** |
| longstring | 0.68 | 0.71 | 0.80 | **0.97** |
| mixed | 0.68 | 0.68 | 0.81 | **0.98** |

### Interpretazione

**Split-half MIN non è un classificatore binario** — è un **estimatore continuo
del grado di corruzione**. La curva di detection sale quasi linearmente dal
60% di flagging a 10% corruzione fino al 99% a 60%+. Utile per ranking.

**Split-half MEAN** è il classificatore binario migliore (+0.007 MCC su VP).

### Raccomandazioni aggiornate

- **Classificazione binaria:** `ReReReRe(data, auto_z=TRUE, variance_penalty=TRUE)` — MCC=0.582
- **Se la corruzione parziale è la preoccupazione principale:** usare
  `score_split_half(data, aggregation="mean")` — MCC=0.589
- **Per score continuo di severity:** `score_split_half(data, aggregation="min")`
  + soglia percentile (es. top-20%)

---

# Aggiornamenti precedenti — 2026-04-02

## Decisione architetturale: due versioni del ReReReRe

A partire dal 1 aprile il codice espone **due funzioni** nel file `ReReReRe.R`:

| Funzione | Metodo | Default | Status |
|----------|--------|---------|--------|
| `ReReReRe()` | **Standard a 2 livelli**: weighted (tutte le coppie, pesate per \|r\|) per ≤60 item, coupled (top-k% per \|r\|) per >60 item | Sì | **Validato** su dati reali |
| `ReReReRe_F()` | **Per-factor**: EFA → coppie within-factor → pesate per \|r\| | No | **Sperimentale** |

### Perché due versioni?

1. **Standard (`ReReReRe`)** vince sui dataset di validazione esterna con ground truth (AUC=0.635 vs EFA-D=0.561, coupled vince 5/6 dataset)
2. **Per-factor (`ReReReRe_F`)** vince su dati simulati (MCC=0.333 vs 0.305 overall) e su questionari corti reali (es. Pennycook)
3. La discrepanza simulazione→reale dipende dalla qualità dell'EFA su dati reali con cross-loadings e strutture complesse

---

## Simulazione ALL VARIANTS (2026-04-02) — Tutte le idee di buone_idee.md testate

**75 condizioni** (5 nF × 3 ipf × 5 reps), N=300, 15% careless, livelli corruzione 10-100%, GT = >50% corruzione. Runtime: 123 minuti.

**13 metodi confrontati:**
- `std` / `std_cf` — standard ReReReRe (con/senza cross-factor baseline)
- `efa_d` / `efa_d_cf` — EFA-D per-factor (con/senza cross-factor baseline)
- `iterative` — EFA iterativa (idea 4): EFA → flag → re-EFA su puliti → ri-score tutti
- `pf_mean` — media dei z per-factor
- `proplow_fz0.5/1.0/1.5/2.0` — % fattori con z < soglia (idea 1, 4 soglie)
- `pf_comb_l0.5/1.0/1.5` — mean_z - λ√var_z (idea 1, 3 lambda)

### Risultati principali: Oracle MCC

| # | Metodo | Oracle MCC |
|---|--------|------------|
| 1 | **iterative** | **0.398** |
| 2 | efa_d | 0.394 |
| 3 | efa_d_cf | 0.390 |
| 4 | std | 0.383 |
| 5 | std_cf | 0.380 |
| 6 | pf_mean | 0.379 |
| 7 | proplow (best) | 0.280 |
| 8 | pf_comb (best) | 0.277 |

### Risultato per fascia di item

| Bin | efa_d | iterative | std | pf_mean | proplow best | Vincitore |
|-----|-------|-----------|-----|---------|-------------|-----------|
| <30 | **0.207** | 0.196 | 0.196 | 0.198 | 0.181 | efa_d |
| 30-60 | **0.281** | 0.279 | 0.261 | 0.258 | 0.199 | efa_d |
| 60-100 | 0.379 | **0.407** | 0.324 | 0.355 | 0.272 | iterative |
| 100-200 | 0.577 | 0.580 | 0.576 | 0.571 | 0.383 | iterative/std (≈pari) |
| >200 | 0.724 | 0.725 | **0.832** | 0.705 | 0.585 | **std** |

### Vincitore per cella (nF × ipf)

| nF\ipf | 3 | 6 | 10 |
|--------|---|---|-----|
| 4 | std (0.185) | efa_d (0.234) | efa_d (0.297) |
| 8 | efa_d_cf (0.224) | pf_mean (0.306) | iterative (0.404) |
| 12 | iterative (0.254) | iterative (0.428) | iterative (0.531) |
| 20 | iterative (0.284) | std (0.542) | iterative (0.698) |
| 30 | iterative (0.390) | efa_d_cf (0.599) | **std (0.832)** |

**L'iterative vince 7/15 celle, specialmente nella fascia 60-200 items.**

### Conclusioni per ogni idea testata

#### Idea 2 (cross-factor baseline): ❌ NESSUN BENEFICIO

| Metodo | Senza CF | Con CF | Differenza |
|--------|----------|--------|------------|
| efa_d | 0.394 | 0.390 | **−0.004** |
| std | 0.383 | 0.380 | **−0.003** |

Il campionamento cross-factor delle coppie random non migliora la performance. Probabilmente perché la quota di coppie accidentalmente within-factor nella baseline è già piccola (~1/nF) e la loro rimozione non compensa il costo di dover fare EFA.

#### Idea 1 — proplow (% fattori con z < soglia): ❌ NETTAMENTE PEGGIORE

MCC massimo 0.280 vs 0.383 per std. Tutte e 4 le soglie testate (0.5, 1.0, 1.5, 2.0) danno risultati simili e inferiori al metodo standard. Il problema è che la binarizzazione (z < soglia → incoerente) perde informazione rispetto allo z-score continuo.

#### Idea 1 — pf_comb (mean_z − λ√var_z): ❌ MOLTO PEGGIORE

MCC 0.061-0.277. La varianza dello z per-factor non aggiunge informazione utile; al contrario, la penalità per alta varianza danneggia la classificazione perché sia i careless che gli attentivi possono avere alta varianza inter-fattore.

#### Idea 1 — pf_mean (media z per-factor): ➖ EQUIVALENTE A STD

MCC 0.379 vs std 0.383. In pratica identico. La media dei z per-factor non è né meglio né peggio del z-score globale.

#### Idea 4 (iterative EFA): ✅ VINCITORE, +4% su std

MCC 0.398 vs std 0.383 (+3.9%). Il vantaggio è concentrato su 60-200 items dove la ri-stima dell'EFA su dati puliti migliora la qualità delle coppie selezionate. A >200 items lo std coupled è meglio (0.832 vs 0.725) perché con tanti items la selezione empirica per |r| è già ottimale.

---

## Raccomandazione aggiornata

| Situazione | Metodo raccomandato |
|------------|---------------------|
| Uso generale / analisi primaria | `ReReReRe()` (standard auto-switch) |
| Questionari corti (≤60 items) | `ReReReRe_F()` (per-factor EFA-D) |
| Questionari medi (60-200 items) | `ReReReRe_F()` con iterative (da implementare) |
| Questionari lunghi (>200 items) | `ReReReRe()` (standard coupled) |

**Nota:** l'iterative EFA attualmente è implementata come `ReReReRe_F_iterative()` separata. Potrebbe essere integrata come opzione in `ReReReRe_F()`.

**Parametri fissi raccomandati (invariati):**
- corProp = 0.03
- z_threshold = 1.5
- auto_z = TRUE per comodità

---

## Metodi EFA testati e scartati come default

### Opzione A — EFA + loading weights (λ_i × λ_j)
54 condizioni simulate. Vince solo a 60-100 item. **Scartata.**

### Opzione D (= EFA-D) — EFA within-factor + |r| weights
360 condizioni simulate. Overall simulazione: EFA-D=0.347 > Std=0.303. Ma su dati reali: EFA-D AUC=0.561 vs Coupled AUC=0.635 (coupled vince 5/6 dataset). **Scartata come default.** Rimane via `ReReReRe_F()`.

---

## Test su Pennycook & Rand 2019

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

**Su questionari corti, EFA-D e weighted migliorano i risultati. Coupled non funziona.**

---

## Nota sulla qualità dei dataset di validazione esterna

I dataset con ground truth di careless responding presentano problemi:
- **Goldammer S1-S3**: 63-67% careless (sperimentale, instructed) — rate irrealisticamente alto
- **Niessen**: ground truth = speed manipulation, non necessariamente inconsistenza
- **Schneider**: ground truth = classe latente (algoritmica, non sperimentale)
- **Schroeders**: il più pulito, ma 40% careless instructed

La differenza tra "vince su simulazione" e "perde su dati reali" potrebbe essere almeno in parte dovuta alla qualità discutibile del ground truth, non solo ai limiti dell'EFA.

---

## File e struttura

### File core (root)

| File | Funzione | Scopo |
|------|----------|-------|
| `ReReReRe.R` | `ReReReRe()` + `ReReReRe_F()` + varianti | Algoritmo principale |
| `Synthetic_Good_Responses_2.R` | `simulated_good_responses()` | Genera dati questionario puliti |
| `Careless_machine_2.R` | `inject_careless()` | Inietta careless responding |

### Simulazioni recenti (root)

| File | Scopo | Status |
|------|-------|--------|
| `Simulation_AllVariants.R` | Confronto 13 metodi (75 cond, 5 reps) | ✅ Completata |
| `Simulation_Realistic_Corruption.R` | Std vs F con corruzione 10-100% | ✅ Completata |
| `sim_variants_results.csv` | Risultati all-variants (12,525 righe) | ✅ |
| `sim_realistic_results.csv` | Risultati corruzione realistica (1,125 righe) | ✅ |

### Report

| File | Contenuto |
|------|-----------|
| `archive/variants_comparison/report.txt` | Report all-variants |
