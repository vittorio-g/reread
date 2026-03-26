# Aggiornamenti sessione 2026-03-26

## Cosa e' stato fatto

### 1. Simulazione di calibrazione nF -> z_threshold
- Script: `Calibration_nF_to_Z.R`
- Variato nF da 2 a 40 (step 1), z_threshold da 0.1 a 3.0 (step 0.2)
- 30 repliche per ogni nF, items/factor fisso a 6, n=300, 10% careless
- Totale: 1,170 chiamate ReReReRe, completate in ~47 minuti
- Risultati in `archive/calibration_nF_z/`

**Risultato chiave:** lo z_threshold ottimale segue una curva a U in funzione di nF:
- nF 2-7: z ~1.0 (ma MCC bassissimo, zona morta)
- nF 8-14: z scende a ~0.5-0.7 (segnale emerge)
- nF 15-20: z stabile ~0.5-0.7 (sweet spot)
- nF 21-28: z risale a ~0.7-1.0
- nF 29-40: z sale a ~1.3-2.1 (segnale forte)

**Grafico:** `archive/calibration_nF_z/plot_calibration_scatter.png` -- scatter cloud con
punti individuali (dimensione = MCC) e media evidenziata in rosso. Molto efficace visivamente.

### 2. Auto-calibrazione implementata in ReReReRe.R
- Nuovo parametro `auto_z=TRUE`: stima automaticamente nF via parallel analysis (`psych::fa.parallel`)
- Lookup del z_threshold ottimale via curva LOESS (span=0.4) fittata sui 39 punti di calibrazione
- La curva LOESS e' embedded direttamente nel codice (nessuna dipendenza da CSV esterni)
- Funzione `get_calibrated_z(nf)` esportata per uso diretto
- Clamping a nF=[2,40] per input fuori range
- Fallback a z=1.5 se la parallel analysis fallisce
- Output arricchito: `z_threshold_used`, `nFactors_detected`
- `flagged` ora usa z_score (non piu' percentile)
- Testato end-to-end su dati simulati: 10 fattori, PA rileva 7, z calibrato a 0.95

### 3. Script di validazione esterna preparato
- `Test_Auto_Z_Validation.R`: confronta auto_z vs z fisso (1.5, 2.0) vs Mahalanobis
- Pronto per i 3 dataset reali (Schneider, Schroeders, Niessen)
- **Non ancora eseguito**: i dataset esterni non sono presenti nel repo dopo il reset

### 4. Organizzazione file
- Risultati spostati in `archive/` (calibrazione, multiverse precedente, benchmark, checkpoint)
- 15 cartelle dataset create in `Dataset/`
- `CLAUDE.md` aggiornato con tutti i risultati e le note di design

## Prossimi passi
- Recuperare i dataset di validazione esterna (Schneider `carersp.csv`, Schroeders `data_mod_resp.csv`, Niessen `Raw_data.sav`)
- Eseguire `Test_Auto_Z_Validation.R` sui dati reali
- Valutare se la parallel analysis sottostima nF in modo sistematico e se serve un correttivo
- Decidere se auto_z diventa il default o resta opt-in
