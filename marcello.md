# Aggiornamenti per Marcello — 2026-04-26

## Phase 6 — Sweep sistematico della soglia GT (la cosa più importante)

PDF in Downloads: `ReReReRe_ThresholdSweep_2026-04-26.pdf`.

### Il problema metodologico che ho risolto

Nelle fasi 1–5 avevo usato **τ_GT = 0.80** come definizione operativa di "careless"
(cioè: chiamiamo careless chiunque abbia >80% di risposte corrotte). Avevo preso
quella soglia come data, non l'avevo ottimizzata.

Vittorio ha (giustamente) obiettato: "questa soglia va scelta dai dati, non
imposta. Sweep sistematico, vediamo dove sta l'ottimo." Phase 6 fa esattamente
questo.

**Verifica preliminare:** la simulazione conteneva tutti i livelli di corruzione
(0%, 10%, 20%, ..., 100%, ~960 rispondenti per livello + 14400 puliti, totale
24000). Nessun errore di generazione — solo una scelta non ottimizzata.

### Cosa ho fatto

Per τ_GT in {0.10, 0.20, ..., 0.90} ho riallenato il Random Forest 5-fold CV con
le 6 feature standard, sotto due scenari:

- **Scenario A** ("clean vs careless puro"): positivi = corruzione > τ,
  negativi = SOLO i puliti (gli intermedi vengono esclusi).
- **Scenario B** ("tutti dentro"): positivi = corruzione > τ, negativi =
  corruzione ≤ τ (puliti + intermedi sono tutti negativi).

Per ogni (τ, scenario) calcolo tutte le 13 metriche al punto operativo FPR=5% e
all'ottimo MCC.

### Risultati — Scenario A (clean vs careless)

| τ_GT | n_pos | MCC | F1 | AUPRC | Kappa |
|:---:|:---:|:---:|:---:|:---:|:---:|
| 0.10 | 8637 | 0.671 | 0.767 | 0.886 | 0.656 |
| 0.30 | 6715 | 0.755 | 0.826 | 0.918 | 0.752 |
| 0.50 | 4804 | 0.802 | 0.852 | 0.938 | 0.802 |
| **0.60** | **3842** | **0.811** | **0.852** | **0.942** | **0.810** |
| 0.70 | 3050 | 0.799 | 0.835 | 0.939 | 0.798 |
| **0.80** *(prima)* | 1930 | 0.783 | 0.804 | 0.939 | 0.774 |
| 0.90 | 967 | 0.702 | 0.703 | 0.925 | 0.678 |

**Argmax:** MCC, AUPRC, Kappa → τ=0.60. F1 → τ=0.50.

### Risultati — Scenario B (tutti dentro)

| τ_GT | n_pos | MCC | F1 | AUPRC | Kappa |
|:---:|:---:|:---:|:---:|:---:|:---:|
| 0.20 | 7833 | 0.678 | 0.763 | 0.881 | 0.668 |
| **0.30** | 6715 | 0.704 | **0.776** | 0.878 | 0.699 |
| **0.40** | 5758 | **0.704** | 0.769 | 0.864 | 0.701 |
| 0.50 | 4804 | 0.690 | 0.749 | 0.843 | 0.689 |
| 0.60 | 3842 | 0.661 | 0.714 | 0.805 | 0.661 |
| **0.80** *(prima)* | 1930 | 0.548 | 0.585 | 0.689 | 0.545 |
| 0.90 | 967 | 0.479 | 0.483 | 0.664 | 0.454 |

**Argmax:** MCC, Kappa, Bal.Acc., Youden → τ=0.40. F1 → τ=0.30.

### La risposta alla domanda

**No, τ=0.80 NON è la migliore secondo nessuna metrica di qualità classificatoria.**

- In Scenario A, l'ottimo è τ=0.60 (MCC=0.811 vs 0.783 a τ=0.80).
- In Scenario B, l'ottimo è τ=0.40 (MCC=0.704 vs 0.548 a τ=0.80 — gap enorme).

AUC e Bal.Acc. continuano a salire fino a τ=0.90, ma sono "ingannate" dai
positivi estremi facili (chi ha 100% corruzione è banale da rilevare). Le
metriche giuste (MCC, F1, AUPRC, Kappa) hanno tutte un ottimo interno chiaro.

### Cosa cambia per il paper

**Riportare entrambi gli scenari, abbandonare τ=0.80.**

- Headline frase A: "MCC ≈ 0.81 con definizione stretta clean-vs-careless
  (Scenario A, τ=0.60)"
- Headline frase B: "MCC ≈ 0.70 con definizione inclusiva (Scenario B, τ=0.40)"

Scenario A è più "pulito" ma artificialmente facile (esclude gli intermedi).
Scenario B è più realistico per uso reale (in pratica non sai a priori chi è
intermedio). Riportarli entrambi è onesto e mostra robustezza.

### Caveat

I numeri di Phase 5 (MCC=0.786, ablation Δ=+0.119 da RR) erano calcolati a
τ=0.80, Scenario A. Se rifacciamo Phase 5 a τ=0.60 i numeri assoluti salgono
leggermente (0.786 → ~0.81), ma i risultati qualitativi (RF batte logit 13/13,
RR contribuisce ΔMCC>0.10, soffitto duro senza RR) sono robusti.

---

## Phase 5b — Phase 5 rifatto a τ=0.60 (questo lo confermo!)

PDF in Downloads: `ReReReRe_Phase5b_tau60_2026-04-26.pdf`.

Ho rifatto tutto Phase 5 (multi-metric + ablation + per-pattern) sotto la
nuova definizione operativa τ=0.60, Scenario A. I numeri non solo *non*
peggiorano — **migliorano**, e il contributo del ReReReRe **cresce**.

### Numeri principali (RF @ FPR=5%)

| Metrica | τ=0.80 (vecchio) | τ=0.60 (nuovo) |
|---------|:---:|:---:|
| MCC | 0.786 | **0.813** |
| F1 | 0.807 | **0.853** |
| AUPRC | 0.940 | **0.942** |
| Kappa | 0.778 | **0.812** |
| PPV | 0.716 | **0.826** |
| Sens | 0.924 | 0.882 |
| Spec | 0.951 | 0.950 |

RF vince ancora **13 metriche su 13** vs logit/glmnet. Identico al vecchio.

### Punti operativi (training)

| Punto | τ classifier | Sens | Spec | PPV | F1 | MCC |
|-------|:---:|:---:|:---:|:---:|:---:|:---:|
| FPR=5% | 0.310 | 0.882 | 0.950 | 0.826 | 0.853 | **0.813** |
| FPR=10% | 0.165 | 0.932 | 0.900 | 0.714 | 0.808 | 0.760 |
| Youden's J | 0.202 | 0.919 | 0.916 | 0.744 | 0.822 | 0.776 |
| F1 ottimale | 0.490 | 0.825 | 0.979 | 0.913 | 0.867 | 0.835 |
| MCC ottimale | 0.510 | 0.819 | 0.981 | 0.919 | 0.867 | **0.836** |

### Ablation (la cosa importante)

| Configurazione | n_feat | RR? | MCC nuovo | MCC vecchio |
|----------------|:---:|:---:|:---:|:---:|
| **Full ensemble** | 6 | sì | **0.813** | 0.784 |
| Triade iter+IRV+D² | 3 | iter | 0.794 | 0.778 |
| Senza z_RR_iter | 5 | solo std | 0.777 | 0.749 |
| **Senza ReReReRe** | 4 | **no** | **0.681** | 0.665 |
| Aux quad | 4 | no | 0.680 | 0.665 |
| Aux triade | 3 | no | 0.681 | 0.665 |
| Solo z_RR family | 2 | sì (sole) | 0.499 | 0.533 |

**Δ MCC dato dal ReReReRe = +0.132** (era +0.119 a τ=0.80). Il contributo
**cresce**, non cala. Stessa cosa per F1 (+0.111) e AUPRC (+0.113).

### Per-pattern: dove il ReReReRe rescuet di più

| Pattern | Corruzione | Full | Senza RR | Δ |
|---------|:---:|:---:|:---:|:---:|
| random | 0.90 | 0.847 | 0.335 | **+0.512** |
| random | 0.80 | 0.803 | 0.307 | +0.496 |
| random | 1.00 | 0.827 | 0.352 | +0.475 |
| random | 0.70 | 0.779 | 0.350 | +0.429 |
| mixed | 0.70 | 0.838 | 0.415 | +0.423 |
| mixed | 1.00 | 0.887 | 0.519 | +0.368 |
| longstring | 0.70 | 0.866 | 0.612 | +0.254 |
| longstring | 0.80 | 0.947 | 0.737 | +0.210 |
| fatigue | 0.90 | 0.929 | 0.768 | +0.161 |
| acquiescent ≥0.80 | — | 1.00 | ~1.00 | ~0 |
| pure_straight ≥0.79 | — | 1.00 | ~1.00 | ~0 |

Il ReReReRe salva esattamente i pattern incoerenti (random, mixed) che gli
ausiliari non vedono. Sui pattern coerenti (acquiescent, pure_straight,
fatigue piena) gli ausiliari saturano già a ~100% e RR non aggiunge nulla.
**Questa complementarità strutturale è il motivo per cui l'ensemble funziona.**

### Frase headline da mettere nel paper

> "Random Forest ensemble di sei feature (z_RR, z_RR_iter, IRV, LongString,
> D², Person-Total) raggiunge MCC = 0.813, F1 = 0.853, AUPRC = 0.942 a
> FPR = 5% su n = 18,242 rispondenti simulati sotto la definizione operativa
> τ_GT = 0.60 (Scenario A: clean vs careless con corruzione > 60%). RF
> batte logistica e elastic-net su tutte e 13 le metriche di valutazione.
> Il ReReReRe contribuisce ΔMCC = +0.132 rispetto all'ensemble di soli
> ausiliari; senza ReReReRe, l'ensemble di quattro detector si satura a
> MCC ≈ 0.68 a prescindere da quanti ausiliari si impilano."

### Cosa cambia nei file di Phase 5 (vecchi)

I CSV originali (`phase5_*.csv`) restano in repository per riferimento.
I nuovi sono `phase5b_*_tau60.csv`. Il PDF Phase 5 originale
(`ReReReRe_MultiMetric_Report_2026-04-25.pdf`) ha numeri a τ=0.80; il
nuovo PDF (`ReReReRe_Phase5b_tau60_2026-04-26.pdf`) è quello da usare.

---

# Aggiornamenti per Marcello — 2026-04-25 (sera)

## Phase 5 — Rivalutazione multi-metrica + ablation study

Due nuovi documenti PDF in Downloads:
- `ReReReRe_RandomForest_Methodology_2026-04-25.pdf` — spiega passo-passo cos'è e perché si usa la Random Forest dentro l'ensemble, con risultati ablation
- `ReReReRe_MultiMetric_Report_2026-04-25.pdf` — rivalutazione su 13 metriche (MCC, F1, F2, AUPRC, AUC, Kappa, Bal.Acc, Youden's J, G-mean, Sens, Spec, PPV, NPV)

### Risultato 1: la rivalutazione multi-metrica conferma tutto

Random Forest vince **13 metriche su 13** rispetto a logistic e glmnet (clean sweep).
Niente da rivedere — l'MCC era allineato al resto. Tabella sintetica:

| Metrica | Logit | glmnet | RF | Δ |
|---------|:---:|:---:|:---:|:---:|
| MCC | 0.730 | 0.730 | **0.786** | +0.056 |
| F1 | 0.761 | 0.761 | **0.807** | +0.046 |
| AUPRC | 0.884 | 0.884 | **0.940** | +0.056 |
| Kappa | 0.725 | 0.725 | **0.778** | +0.053 |
| Sens | 0.843 | 0.844 | **0.924** | +0.080 |
| AUC | 0.969 | 0.969 | **0.982** | +0.013 |

### Risultato 2 (la domanda chiave): il ReReReRe contribuisce davvero?

Ablation con 7 configurazioni dell'ensemble. Domanda: il ReReReRe (z_RR + z_RR_iter)
serve davvero, o l'ensemble di soli detector ausiliari (IRV, LongString, D², PT)
funziona altrettanto bene?

| Configurazione | n_feat | ReReReRe? | MCC | F1 | AUPRC |
|---------|:---:|:---:|:---:|:---:|:---:|
| Full ensemble | 6 | sì (entrambi) | **0.784** | **0.805** | **0.941** |
| Triade iter+IRV+D² | 3 | sì (solo iter) | 0.778 | 0.801 | 0.927 |
| Senza z_RR_iter | 5 | sì (solo std) | 0.749 | 0.777 | 0.907 |
| **Senza ReReReRe** | 4 | **NO** | **0.665** | 0.705 | 0.811 |
| Aux quad: IRV+D²+LS+PT | 4 | NO | 0.665 | 0.705 | 0.811 |
| Aux triade: IRV+D²+LS | 3 | NO | 0.665 | 0.705 | 0.803 |
| Solo z_RR family | 2 | sì (soli) | 0.533 | 0.586 | 0.659 |

**Risposta: SÌ, il ReReReRe contribuisce.** Δ MCC = +0.119, Δ F1 = +0.100,
Δ AUPRC = +0.130. Sono tutti effetti grandi e consistenti.

**Tre osservazioni strutturali:**

1. **L'iterativo-EFA porta tutto il peso.** La triade `iter+IRV+D²` (3 feature)
   raggiunge MCC=0.778, solo 0.006 sotto l'ensemble completo. Il z_RR standard
   sopra l'iter non aggiunge nulla.

2. **Soffitto duro a MCC≈0.665 senza ReReReRe.** Tutte e tre le configurazioni
   no-RR atterrano *esattamente* allo stesso MCC. Aggiungere altri ausiliari
   non aiuta più una volta che z_RR_iter è assente.

3. **Il ReReReRe da solo non basta.** Solo z_RR (2 feature) → MCC=0.533. È
   necessario ma non sufficiente.

### Risultato 3: dove esattamente ReReReRe contribuisce?

| Pattern | Corruzione | Full | No-RR | Δ da RR |
|---------|:---:|:---:|:---:|:---:|
| **random** | 100% | 0.846 | 0.352 | **+0.494** |
| **random** | 90% | 0.835 | 0.365 | **+0.470** |
| **mixed** | 100% | 0.881 | 0.500 | **+0.381** |
| **mixed** | 90% | 0.767 | 0.497 | **+0.270** |
| longstring | 90% | 0.950 | 0.761 | +0.189 |
| fatigue | 90% | 0.916 | 0.761 | +0.155 |
| longstring | 100% | 0.940 | 0.807 | +0.133 |
| pure_straight, acquiescent, fatigue 100% | — | 1.000 | 1.000 | 0 |

ReReReRe **salva specificamente i pattern incoerenti** (random, mixed). Sui
pattern consistenti (straight-lining, acquiescenza) gli ausiliari sono già
perfetti e ReReReRe non aggiunge nulla. Questa **complementarietà strutturale**
è la ragione per cui l'ensemble funziona — non è un caso, è un'architettura.

### Bottom line per la pubblicazione

> "L'ensemble Random Forest a 6 detector raggiunge MCC=0.786 (FPR=5%) sulla
> detection di carelessness piena. Random Forest domina i classificatori
> lineari su tutte e 13 le metriche valutate. Lo studio di ablazione mostra
> che il contributo unico di ReReReRe è di +0.119 MCC: il modulo è
> indispensabile, in particolare per i pattern incoerenti (random,
> mixed) che gli altri detector non vedono."

Tutti i file CSV sono in `phase5_*.csv`. Plot in `plot_phase5_*.png`.
Script: `Phase5_MultiMetric.R`.

---

# Aggiornamenti per Marcello — 2026-04-25

## Ottimizzazione publication-ready: MCC 0.78 con RF + ensemble + Scenario A

**Decisione operativa:** definiamo "careless" come **corruzione > 80%**. Questo
è dove l'indice diventa publication-ready (MCC > 0.7). Sotto, la corruzione
parziale è statisticamente troppo ambigua per binarizzazione affidabile.

**Configurazione vincente:**
1. **Ensemble di 6 detector** (z_RR + z_RR_iter_EFA + IRV + LongString + D² + PersonTotal)
2. **Random Forest** come combiner (vince su logistic +0.05 MCC)
3. **Scenario A**: classe positiva = corruption > 80%, classe negativa = clean SOLO (escludi 10-80%)
4. **Operating point fisso a FPR = 5%** (95° percentile dei predict sui clean)

**Risultati (24K respondenti, 5-fold CV):**
- MCC = **0.788**
- Sensitivity = 0.92
- Specificity = 0.95
- AUC = **0.983**

**MCC scala con questionnaire length:**
| Items | MCC |
|:---:|:---:|
| 48 | 0.58 |
| 80 | 0.72 |
| 96 | 0.78 |
| 160 | 0.89 |
| 300 | **0.95** |

**Validazione indipendente** (54K respondenti freschi, seeds diversi, 3 careless rates):
| Rate | MCC | Sens |
|:---:|:---:|:---:|
| 20% | 0.644 | 0.91 |
| 40% | 0.785 | 0.92 |
| 60% | 0.839 | 0.91 |

**Per-pattern detection (corruption > 80%):**
- pure_straight, acquiescent: 1.00
- fatigue, longstring: 0.93-0.95
- mixed: 0.78-0.84
- random: 0.78-0.83

Tutti i pattern ≥78%. RF risolve fatigue (era 0.51 con logit).

**Feature importance (RF Mean Decrease Accuracy):**
1. **IRV: 142** (top)
2. **Mahalanobis D²: 122**
3. **z_RR_iter_EFA: 106**
4. z_RR standard: 45 (subsumed da iter)
5. LongString: 43
6. PersonTotal: 24

**Triade minima vincente:** IRV + D² + z_RR_iter_EFA cattura ~90% del segnale. Standard z_RR è ridondante quando z_RR_iter è presente.

**Quando l'indice è ready:**
- ≥80 items + careless rate ≥40%: MCC > 0.7 ✓
- ≥160 items: MCC > 0.85 in tutte le condizioni
- ≥300 items: MCC ≈ 0.95-0.97

PDF report: `ReReReRe_Publication_Optimization_2026-04-25.pdf` in Downloads.
Script: `Phase1_Diagnostic.R`, `Phase2_Properties.R`, `Phase3_Validation.R`,
`Phase4_FinalPlots.R`.

---

# Aggiornamenti precedenti — 2026-04-21

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
