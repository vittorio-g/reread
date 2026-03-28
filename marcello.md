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
