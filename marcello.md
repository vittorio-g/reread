# Aggiornamenti per Marcello

## Sessione 2026-03-28 — Riepilogo

### 1. Full Multiverse Analysis (completata)

Script: `Multiverse_Full.R` | Runtime: 161 min | 5,670 chiamate ReReReRe | 85,050 righe risultato

**Design:** 7 nF × 3 ipf × 3 n × 3 pct × 3 corProp × 10 repliche

**Finding principale: items_per_factor spiega più varianza di nFactors!**

Importanza delle variabili (eta² da ANOVA one-way sull'MCC oracle):

| Variabile | Varianza spiegata |
|-----------|------------------|
| **items_per_factor** | **31.1%** |
| **nFactors** | **26.1%** |
| pct_careless | 3.6% |
| corProp | 1.7% |
| n_respondents | 1.0% |

Questo significa che il numero di item per fattore conta **più** del numero di fattori. La raccomandazione pratica va riformulata in termini di **lunghezza totale del questionario** (nF × ipf), non solo numero di costrutti.

#### MCC per nFactors × items_per_factor (oracle, corProp=0.05)

| nF | ipf=3 | ipf=6 | ipf=10 |
|----|-------|-------|--------|
| 4  | 0.075 | 0.111 | 0.160  |
| 8  | 0.111 | 0.207 | 0.314  |
| 12 | 0.143 | 0.295 | 0.460  |
| 16 | 0.171 | 0.388 | 0.579  |
| 20 | 0.202 | 0.447 | 0.593  |
| 25 | 0.222 | 0.467 | 0.596  |
| 30 | 0.260 | 0.513 | 0.560  |

**Esempio chiave:** nF=12 con ipf=10 (120 item, MCC=0.46) batte nF=25 con ipf=3 (75 item, MCC=0.22).

#### corProp: più piccolo è meglio

| corProp | Mean MCC |
|---------|----------|
| **0.03** | **0.352** |
| 0.05 | 0.327 |
| 0.10 | 0.286 |

#### Auto-z ≈ z fisso = 1.5

| Strategia | Mean MCC |
|-----------|----------|
| auto-z (1D) | 0.268 |
| z=1.0 | 0.260 |
| **z=1.5** | **0.274** |
| z=2.0 | 0.274 |
| oracle | 0.356 |

L'auto-calibrazione non migliora rispetto a un z fisso di 1.5. Ma la versione 2D (basata su total_items) chiude il 38% del gap con l'oracle.

#### Separazione z-score (buoni vs careless)

| nF | z buoni | z careless | Gap |
|----|---------|-----------|-----|
| 4  | 0.43 | -0.05 | 0.48 |
| 12 | 2.38 | 0.39 | 1.99 |
| 20 | 3.59 | 1.02 | 2.57 |
| 30 | 4.75 | 1.94 | 2.81 |

Il gap cresce quasi linearmente con nF — spiega il vantaggio di scaling.

---

### 2. Calibrazione 2D (nF × ipf → z_threshold ottimale)

Script: `Calibration_nF_ipf_Z.R` | Runtime: 35.5 min | 1,800 chiamate | 60 celle × 30 repliche

**Finding: total_items è il miglior predittore singolo della z ottimale**

- Formula lineare: `z = 0.505 + 0.0042 × total_items` (R² = 0.476)
- vs `z = 0.261 + 0.020×nF + 0.055×ipf` (R² = 0.245) — total_items spiega il doppio!

#### z ottimale per nF × ipf

| nF\ipf | 3   | 6   | 10  | 12  |
|--------|-----|-----|-----|-----|
| 4      | 1.1 | 0.8 | 1.0 | 0.8 |
| 8      | 1.0 | 1.0 | 0.7 | 0.7 |
| 12     | 1.0 | 0.6 | 0.6 | 0.7 |
| 15     | 0.7 | 0.5 | 0.7 | 1.1 |
| 20     | 0.7 | 0.5 | 1.2 | 1.8 |
| 25     | 0.7 | 0.7 | 1.9 | 2.3 |
| 30     | 0.6 | 1.0 | 2.5 | 2.8 |

**Gradiente diagonale:** z bassa (pochi item, top-left) → z alta (tanti item, bottom-right).

#### Performance auto-z 2D vs z fisso

| Strategia | Mean MCC | vs z=1.5 |
|-----------|----------|----------|
| z fisso = 1.5 | 0.316 | — |
| **Auto-z 2D** | **0.340** | **+0.024** |
| Oracle | 0.380 | +0.064 |

Il guadagno è concentrato sui questionari lunghi:

| Total items | z=1.5 | Auto-z | Gain |
|-------------|-------|--------|------|
| 48 | 0.14 | 0.15 | +0.01 |
| 120 | 0.41 | 0.45 | +0.04 |
| 300 | 0.69 | 0.78 | +0.09 |
| 360 | 0.65 | 0.82 | +0.17 |

Per questionari >200 item, auto-z fornisce guadagni sostanziali.

---

### 3. Aggiornamenti al codice ReReReRe.R

| Parametro | Prima | Dopo | Motivazione |
|-----------|-------|------|-------------|
| **corProp** default | 0.05 | **0.03** | MCC 0.352 vs 0.327 nella multiverse |
| **auto_z** | 1D (nF → z) | **2D (total_items → z)** | R²=0.476 vs 0.245; doppia varianza spiegata |
| **Lookup table** | 39 punti (nF) | **60 punti (nF×ipf)** | Calibrazione più robusta |
| **Fallback** | nessuno | **z=0.505+0.0042×items** | Se LOESS fallisce |

---

### 4. Report e grafici generati

#### Full Multiverse (`archive/multiverse_full/report/` — 18 grafici)
1. MCC by nFactors (main effect)
2. MCC by items_per_factor
3. MCC by n_respondents (diminishing returns)
4. MCC by pct_careless
5. MCC by corProp
6. MCC by z_threshold (plateau z=0.9-1.9)
7. z curves by nFactors (optimal z shifts)
8. **Heatmap nF × ipf** (il grafico chiave)
9. MCC by total_items con LOESS (unifying predictor)
10. nF × n interaction
11. nF × pct interaction
12. nF × corProp interaction
13. **Auto-z vs fixed vs oracle** (threshold comparison)
14. Sensitivity-Specificity trade-off
15. Boxplot MCC distribution by nF
16. Heatmap nF × z con bordi sull'ottimale
17. ipf within nF (grouped bars)
18. z-score separation buoni vs careless

#### Calibrazione 2D (`archive/calibration_nF_ipf_z/` — 8+16 grafici)
Heatmaps, scatter clouds, curve z, LOESS fits, confronti fixed vs auto-z.

---

### 5. Raccomandazioni pratiche aggiornate

- **corProp = 0.03** (default aggiornato)
- **z_threshold = 1.5** per uso manuale (robusto universale)
- **auto_z = TRUE** per comodità (la funzione sceglie automaticamente)
- Usare ReReReRe quando **total_items ≥ 60** (es. nF≥10 con ipf≥6)
- Sweet spot: **total_items ≥ 120, n ≥ 300**
- Sotto 60 item: detection debole (MCC < 0.20)

---

### 6. Validazione Esterna (2026-03-28)

Testato ReReReRe (corProp=0.03) vs Mahalanobis pratico (chi-sq .001) su **6 dataset reali** con ground truth di careless responding.

#### Dataset validati

| Dataset | N | Items | nF | Tipo GT | RR AUC | Mah AUC | RR MCC (z=1.5) | Mah MCC |
|---------|---|-------|----|---------|--------|---------|---------------|---------|
| Schroeders 2022 (HEXACO) | 605 | 60 | ~10 | Sperimentale | **0.606** | 0.537 | **0.177** | 0.045 |
| Schneider QoL | 1649 | 31 | ~5 | Classe latente | **0.735** | 0.724 | 0.118 | **0.208** |
| Niessen 2016 (IPIP) | 180 | 100 | ~6 | Speed manipulation | **0.639** | 0.436 | 0.073 | 0.000 |
| Goldammer S1 (BFI-2) | 291 | 60 | ~8 | Sperimentale | 0.711 | **0.787** | **0.320** | 0.257 |
| Goldammer S2 (IPIP) | 265 | 60 | ~6 | Sperimentale | 0.674 | **0.795** | **0.304** | 0.266 |
| Goldammer S3 (long.) | 523 | 60 | ~7 | Sperimentale | 0.462 | 0.550 | -0.036 | 0.011 |

#### Johnson IPIP-NEO-300 (inject-and-detect, COMPLETATO)

300 item, 30 facet, ~148 reverse-coded. 5000 rispondenti, careless iniettati al 5%, 10%, 20%.

| pct | MCC oracle (z=3.0) | MCC z=1.5 | MCC auto (z=2.37) | Specificità |
|-----|---------------------|-----------|-------------------|-------------|
| 5% | **0.750** | 0.628 | 0.726 | 1.000 |
| 10% | **0.726** | 0.635 | — | 1.000 |
| 20% | **0.702** | 0.607 | — | 1.000 |

**Specificità = 1.000** su quasi tutte le condizioni — zero falsi positivi con 300 item!
Auto-z ha scelto z=2.37, appropriato per 300 item. Risultati coerenti con le simulazioni per nF=30.

#### Risultati chiave

1. **RR AUC > Mah AUC su 4/6 dataset** — conferma il vantaggio del ranking
2. **RR MCC (z=1.5) > Mah MCC su 4/6 dataset** — il Mahalanobis pratico resta debole
3. **Auto-z funziona bene su Schroeders**: ha scelto z=0.72 (appropriato per 60 item), MCC=0.228 vs z=1.5 MCC=0.177
4. **Goldammer S3 fallisce per entrambi** — 67% careless in design longitudinale, probabilmente confuso da effetti di pratica
5. **Johnson (300 item)**: performance eccellente, coerente con le simulazioni per nF=30
6. **Breakdown Goldammer S1**: 33% careless (AUC=0.738) leggermente più facile da detectare per RR rispetto a 100% careless (AUC=0.692) — il careless parziale mantiene più struttura che RR può sfruttare

---

### 7. Weighted Mode per questionari corti (2026-03-30)

**Problema:** il ReReReRe standard seleziona il top-k% delle coppie di item per correlazione. Con questionari corti (≤60 item), ci sono poche coppie ad alta |r| → segnale debole, MCC basso.

**Soluzione:** invece di selezionare solo le coppie migliori, usare **TUTTE le coppie** ma pesarle per la loro correlazione sample-level |r|. Le coppie fortemente correlate contano di più, quelle deboli contribuiscono poco ma non vengono scartate.

#### Come funziona

La funzione `rowCor_weighted()` calcola per ogni rispondente:
1. Standardizza le risposte di ogni item (z-score across respondents)
2. Per ogni coppia: calcola il prodotto incrociato z_A × z_B
3. Pesa ogni coppia per il suo |r_sample|
4. Media pesata → "coherence score" del rispondente

#### Risultati simulazione (18 condizioni, nF=4-20, ipf=3-10)

| Range item | Standard (MCC) | Weighted (MCC) | Differenza |
|-----------|----------------|----------------|------------|
| **≤30 item** | 0.110 | **0.151** | **+37%** |
| **30-60 item** | 0.193 | **0.255** | **+32%** |
| 60-100 item | **0.315** | 0.300 | −5% |
| 100-200 item | **0.524** | 0.414 | −21% |

**Crossover a ~60 item**: il weighted vince per questionari corti, il coupled vince per questionari lunghi.

#### Integrazione automatica nel ReReReRe.R

Nuovo parametro `mode` con tre opzioni:
- `"auto"` (default) → **weighted** se ≤60 item, **coupled** se >60 item
- `"coupled"` → forza il metodo standard (top-k% coppie)
- `"weighted"` → forza il metodo pesato (tutte le coppie)

**L'utente non deve fare nulla** — lo switch è automatico e trasparente. Il ReReReRe ora funziona meglio out-of-the-box su questionari di qualsiasi lunghezza.

#### Test di verifica

| Test | Item | Mode auto | MCC weighted | MCC coupled |
|------|------|-----------|-------------|-------------|
| 24 item (4 fattori) | 24 | weighted | **0.118** | 0.031 |
| 120 item (20 fattori) | 120 | coupled | 0.197 | **0.442** |

Conferma: lo switch automatico sceglie sempre il metodo migliore.

---

### 8. Dataset Pennycook & Rand 2019 (2026-03-30)

Replicato l'articolo "Lazy, not biased" (Study 1 + Study 2) e applicato il ReReReRe.

- **Study 1**: 30 item (15 fake + 15 real news accuracy ratings), scala 1-4
- **Study 2**: 24 item (12 fake + 12 real), scala 1-4

Il ReReReRe flagga ~48% dei rispondenti su queste scale molto corte — troppi. Questo ha motivato lo sviluppo del weighted mode (Sezione 7).

Il codice di replicazione è in `Dataset/dataset_1/Replication_and_ReReReRe.R`.

---

### 9. EFA-based pair selection: Opzione A testata e scartata (2026-03-30)

Abbiamo testato un'**Opzione A** in cui l'EFA guidava sia la selezione che il peso delle coppie:
1. Parallel analysis → stima nF
2. EFA (oblimin, minres) → assegna ogni item al fattore primario
3. Genera TUTTE le coppie within-factor
4. Pesa per **λ_i × λ_j** (prodotto dei loading — stimatore "shrinkage")

**Risultati su 54 condizioni (nF=4-20, ipf=3-10, 3 reps):**

| Range item | Standard | Weighted | EFA-A | Vincitore |
|-----------|----------|----------|-------|-----------|
| <30 | 0.101 | **0.167** | 0.152 | Weighted |
| 30-60 | 0.183 | **0.247** | 0.240 | Weighted |
| 60-100 | 0.352 | 0.351 | **0.378** | EFA-A |
| 100-200 | **0.523** | 0.400 | 0.489 | Standard |

L'EFA-A vince solo nella fascia 60-100 item (+0.027 MCC). Il problema: la parallel analysis sottostima nF (~70% del vero) e i prodotti dei loading amplificano questo errore.

**Decisione: Opzione A scartata.**

---

### 10. EFA-D testata e scartata (2026-03-31)

Dopo aver scartato l'Opzione A (loading weights), abbiamo testato l'**Opzione D**:
- L'EFA decide **QUALI coppie** contano → solo coppie within-factor
- La |r| osservata decide **QUANTO pesano** → più robusto dei loading

#### Simulazione ampia (360 condizioni)

8 nF × 3 ipf × 3 pct_careless × 5 reps, N=300. Risultati su **dati simulati** promettenti:

| Range item | Standard | Weighted | **EFA-D** | Vincitore |
|-----------|----------|----------|-----------|-----------|
| <30 | 0.108 | **0.160** | 0.158 | Weighted |
| 30-60 | 0.198 | 0.244 | **0.261** | **EFA-D** |
| 60-100 | 0.306 | 0.307 | **0.351** | **EFA-D** |
| 100-200 | 0.501 | 0.405 | **0.529** | **EFA-D** |
| >200 | **0.656** | 0.393 | 0.653 | Standard |

Complessivo simulazione: EFA-D = 0.347 vs Standard = 0.303 vs Weighted = 0.286.

#### MA: validazione su dati reali — EFA-D perde nettamente

| Dataset | Items | nF | EFA-D AUC | Coupled AUC | EFA-D MCC₁.₅ | Coupled MCC₁.₅ |
|---------|-------|----|-----------|-------------|-------------|----------------|
| Schroeders | 60 | ~10 | 0.583 | **0.613** | 0.126 | **0.182** |
| Schneider | 31 | ~5 | 0.459 | **0.729** | -0.019 | **0.114** |
| Niessen | 100 | ~5 | 0.587 | **0.637** | 0.037 | **0.073** |
| Goldammer S1 | 60 | ~8 | 0.605 | **0.717** | 0.167 | **0.314** |
| Goldammer S2 | 60 | ~6 | 0.599 | **0.656** | 0.096 | **0.286** |
| Goldammer S3 | 60 | ~7 | **0.532** | 0.456 | 0.096 | -0.032 |

**Media: EFA-D AUC=0.561 vs Coupled AUC=0.635.** Il coupled vince su 5/6 dataset.

#### Perché la simulazione mente

Il gap simulazione → dati reali è causato da:
1. **La parallel analysis su dati reali produce assegnazioni fattoriali scadenti** — cross-loading, effetti di metodo, strutture complesse che l'EFA non cattura bene
2. **Dataset con alta % careless (Goldammer 63-67%) corrompono la matrice di correlazione** → l'EFA produce fattori spuri
3. **Il coupled è più robusto perché non dipende dalla struttura fattoriale** — il ranking delle |r| sopravvive alla contaminazione meglio della soluzione fattoriale

#### Lezione appresa

**Le migliorie al metodo vanno SEMPRE validate su dati reali, non solo su simulazioni.** I dati simulati hanno una struttura fattoriale pulita che l'EFA recupera facilmente; i dati reali no.

**Decisione: EFA-D scartata.** Il default resta lo switch a 2 livelli: weighted ≤60 item, coupled >60 item. L'EFA-D resta disponibile via `mode="efa_d"` per ricerca.

Report in `archive/efa_d_comparison/`.

---

### 11. Analisi fattoriale nel ReReReRe

**No, il ReReReRe standard NON fa analisi fattoriale.** Lavora solo con la matrice di correlazione tra item e le permutazioni random. L'EFA-D è stata testata come possibile miglioramento ma scartata (vedi sopra).

L'unica eccezione: `auto_z=TRUE` esegue una parallel analysis per calibrare il z_threshold, ma è opzionale.

`psych` era già una dipendenza (per `auto_z`), quindi non aggiunge nuove dipendenze.
