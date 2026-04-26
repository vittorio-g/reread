# ReReReRe — Analisi critica del paper v2 e piano di mitigazione

Data: 2026-04-26 (notte). Autore: Claude. Scopo: enumerare le critiche
metodologiche più probabili da parte di reviewer e implementare
contromisure per quelle rilevanti.

## 1. Catalogo critiche (in ordine di severità)

### S1 — Circolarità "train e test su stessa distribuzione simulata" (alta)
La RF è allenata su 5-fold CV all'interno di ciascuna run, ma ogni run
genera dataset con la stessa distribuzione (lavaan CFA, loadings 0.2-0.8,
6 pattern di careless). La CV verifica generalizzazione dentro un dataset,
non fra dataset. Reviewer ovvio: "state misurando MCC sulle stesse
condizioni viste in training."

**Mitigazione:** holdout cross-run e cross-condizione.
- A1. Cross-rep holdout: per ogni cella (size × rate), train su rep 1-8,
  test su rep 9-12. Pool tra reps.
- A2. Cross-size holdout: train su sizes 30-150, test su 200-300 (e
  viceversa). Verifica trasferibilità tra dimensioni.
- A3. Cross-pattern holdout: train su 4 pattern (random, longstring,
  acquiescent, fatigue), test su 2 (mixed, pure_straight). Verifica che
  la RF non abbia memorizzato firme specifiche.

### S2 — Calibrazione FPR=5% richiede di conoscere chi è clean (alta)
Il punto operativo è il 95° percentile delle predizioni RF sui SOLI
respondenti clean. Nessun ricercatore in produzione ha accesso a questo
oracolo. Spec=0.95 è quindi quasi tautologica.

**Mitigazione:** B1. Calibrazione quantile-based che non assume label.
- Variante 1: threshold = quantile 95° su tutti i respondenti (non solo
  clean). Sotto-stima della specificità ma realistica.
- Variante 2: threshold = quantile (1 - r̂) dove r̂ è una stima del rate
  attraverso il numero di score con z_RR_iter > 2. Self-bootstrapping.
- Variante 3: "fixed-threshold" deployabile: τ=0.5 sulla probabilità RF.
- Confrontare MCC e Δ MCC sotto le tre varianti.

### S3 — Mancanza di intervalli di confidenza (media)
Δ MCC = +0.137 (300 items, Scen A) vs +0.027 (30 items) — ma senza CI è
difficile sostenere la significatività delle differenze.

**Mitigazione:** C1. Bootstrap percentile CI 95% sull'effetto Δ MCC della
ablazione, per ciascuna size. 1000 bootstrap.

### S4 — Definizione di "careless" non univoca (media)
Phase 6 ha sweep su τ_GT e selezionato 0.60 (Scen A) / 0.40 (Scen B). Ma
queste sono ancora scelte definitorie. Reviewer può chiedere: cosa succede
a τ=0.50 o τ=0.30?

**Stato:** già coperto in Phase 6 (PDF "ThresholdSweep_2026-04-26"). Non
serve ri-fare. Aggiungere riferimento a Phase 6 nel testo del v2.

### S5 — Generalizzazione a scale Likert miste (bassa-media)
Tutto generato a 7-point. PISA ha 4/5/6-point misti. Già discusso nel CLAUDE.md
e già documentato come limitazione (raccomandazione: rescaling preprocessing).

**Stato:** documentato come limitazione. Non agire ora.

### S6 — Confronto solo con metodi single-feature (bassa-media)
Confronto interno fra ablazioni dell'ensemble, non con metodi combinati
pubblicati (es. careless::longstring + careless::IRV usati congiuntamente
da Curran 2016, o il PRP di Reise 2003).

**Stato:** out of scope per stanotte. Aggiungere come limitazione.

### S7 — Validazione esterna dell'ensemble non fatta (media-alta)
Su dataset reali (Schroeders, Schneider, Goldammer) abbiamo confrontato
z_RR singolo vs Mahalanobis singolo — non la RF completa con tutti i 6
features. Mancanza importante.

**Mitigazione:** D1. Calcolare le 6 feature standard sui dataset esterni
con GT (Schroeders, Schneider, Goldammer S1/S2). Allenare RF (con CV)
sulla GT reale. Confrontare con z_RR singolo e Mah singolo. Costo medio.

### S8 — Hyperparameter ottimizzazione (bassa)
RF di default (200 alberi, mtry=sqrt). Tuning potrebbe migliorare ma
solitamente di poco. Non agire.

### S9 — Class imbalance estremo non testato (bassa)
Real surveys spesso hanno 1-3% careless, non 10-60%. Nostro
range 10-60% può non rappresentare il caso più difficile.

**Stato:** Phase 6 e v2 hanno il rate=10% come limite inferiore. Aggiungere
nota di limitazione.

### S10 — Solo n=500 per dataset (bassa)
Real surveys variano da n=100 a n=10000. n=500 fisso è una scelta.

**Stato:** Multiverse principale ha già esplorato n=50-1000. Coperto.

## 2. Piano di esperimenti (priorità)

| # | Critica | Esperimento | Priorità | Stima |
|---|---------|-------------|----------|-------|
| A1 | S1 | Cross-rep holdout | ALTA | 5 min |
| A2 | S1 | Cross-size holdout | ALTA | 10 min |
| A3 | S1 | Cross-pattern holdout | ALTA | 10 min |
| B1 | S2 | Calibrazione realistica (3 varianti) | ALTA | 5 min |
| C1 | S3 | Bootstrap CI sull'ablation | ALTA | 10 min |
| D1 | S7 | Validazione esterna RF | MEDIA | 30 min |

Totale stimato: ~70 minuti compresa la rigenerazione del PDF.

Output:
- `robustness_v2_results.csv` — risultati esperimenti A1-A3, B1
- `bootstrap_ci_v2.csv` — risultati C1
- `external_rf_v2.csv` — risultati D1 (se eseguito)
- `ReReReRe_Article_v2.1_2026-04-26.pdf` — PDF aggiornato con sezione
  "Robustness analyses (added in revision)"
