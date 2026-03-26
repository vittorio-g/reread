# External Datasets for ReReReRe Validation

Downloaded 2026-03-25. All datasets contain some form of known/labeled careless respondents.

## Summary Table

| # | Dataset | N | Items | Ground Truth Type | Scales/Domains | Suitability for RR |
|---|---------|---|-------|-------------------|----------------|-------------------|
| 01 | Schroeders et al. 2022 | ~1300 (605 experimental) | ~200+ | **Experimental induction** (instructed careless vs diligent) | Personality, attitudes, misc | **HIGH** — many items, clean labels |
| 02 | Kopitar & Stiglic 2023 | same data as 01 | same | Same (reanalysis with SHAP) | Same | Same data, different analysis |
| 03 | Żółtak et al. 2024 | ~4241 | ~120+ grid items | **Experimental induction** + cursor trajectories | Personality (Big Five) + attitudes | **HIGH** — experimental labels + paradata |
| 04 | Bloy et al. 2025 | ~606 (Study 1) + 1242 (IDs) | 17 gibberish + survey items | **Gibberish scale** (ML-based classification) | Various (gibberish + real items) | **MEDIUM** — novel approach, small N |
| 05 | Brühlmann et al. 2020 | ~394 | 163 (semicolon-delim) | **CrowdFlower quality flags** + self-report + bogus items | PANAS, BFI, needs satisfaction, misc | **HIGH** — rich multi-scale, multiple CR indicators |
| 06 | Schneider et al. (QoL) | ~1987 | ~58 | **Latent class** (careless subgroups) | Depression, pain, fatigue, cognition (PROMIS) | **MEDIUM** — 4 domains, clinical data |
| 07 | Kuang et al. 2025 | ~24,292 | PHQ-9 + PSS + ISI + GAD-7 | **Page time index** (< 2s/item = careless) | Depression, stress, insomnia, anxiety | **LOW for RR** — only 4 short scales (7-14 items each) |
| 08 | Niessen et al. 2016 | unknown (SPSS) | unknown (SPSS) | **Experimental** (speed instructions) | Personality (web questionnaire) | **NEEDS EXPLORATION** — SPSS format |
| — | gjkev (unnamed project) | unknown (SPSS) | unknown (SPSS) | **Longitudinal** (response style changes) | unknown | **NEEDS EXPLORATION** — SPSS format |

## Detailed Notes

### 01_Schroeders_2022 — Detecting CR with Gradient Boosted Trees
- **Paper:** Schroeders, Schmidt & Gnambs (2022). Detecting Careless Responding in Survey Data Using Stochastic Gradient Boosting. *Educational and Psychological Measurement*.
- **Source:** https://osf.io/mct37/
- **Design:** Web-based experiment. Participants randomly assigned to **diligent** vs **careless** responding conditions. Careless group was instructed to respond without reading items carefully.
- **Files:**
  - `careless.csv` — 1302 rows, full dataset with all variables + condition labels
  - `data_mod_resp.csv` — 605 rows, item responses only
  - `data_mod_rt.csv` — response times
  - `data_mod_resp_rt.csv` — responses + times combined
  - `Variables_Listing.pdf` — codebook
- **Ground truth:** `IN01_CP` column indicates experimental condition (careless/diligent).
- **Items:** ~200+ questionnaire items across multiple scales.
- **ReReReRe suitability:** **Excellent.** Large item pool, clean experimental labels, widely cited benchmark.

### 02_Kopitar_2023 — SHAP Analysis of Careless Respondent Characteristics
- **Paper:** Kopitar & Stiglic (2023). Using heterogeneous sources of data and interpretability of prediction models to explain the characteristics of careless respondents in survey data. *Scientific Reports*.
- **Source:** https://github.com/lkopitar/Careless_SHAP
- **Design:** Reanalysis of Schroeders et al. data with SHAP explainability.
- **Files:** Pre-processed train/test splits in .rds format, pre-trained GBM models.
- **Note:** Same underlying data as 01, but with derived features (longstring, response time indices, etc.) already computed. Useful as reference for feature engineering.

### 03_Zoltak_2024 — Cursor Trajectories for Careless Detection
- **Paper:** Żółtak, Pokropek & Muszyński (2024). Identifying Careless Responding in Web-Based Surveys: Exploiting Sequence Data from Cursor Trajectories and Approximate Areas of Interest. *Zeitschrift für Psychologie*.
- **Source:** https://osf.io/sre9g/
- **Design:** Web experiment with **experimentally induced inattentiveness**. Includes item responses, cursor movement data, and response times.
- **Files:**
  - `auxiliary-data-wide.csv` — 4241 respondents, wide format
  - `auxiliary-data-long.csv` — long format (28352 rows)
  - Codebooks for both
  - `results-to-analysis.rds` — preprocessed R data
- **Ground truth:** Experimental condition (attentive vs inattentive). Also includes computed carelessness indicators from cursor data.
- **ReReReRe suitability:** **Very good.** Experimental labels + rich paradata. Need to extract item responses from the auxiliary data.

### 04_Bloy_2025 — Gibberish Scale for Careless Detection
- **Paper:** Bloy, Resheff, Kluger & Malovicki-Yaffe (2025). Identifying Careless Survey Respondents Through Machine Learning Using Responses to a Gibberish Scale. *Advances in Methods and Practices in Psychological Science*.
- **Source:** https://osf.io/32x9k/
- **Files:**
  - `gibberish.csv` — 606 respondents, 17 columns (Q1_1 through Q3_8 + Duration)
  - `Study1_IDs.csv` — 1242 rows, respondent IDs with careless labels
  - `code_Study1.R` — analysis code
  - `GIBML_ForUse_7Point.R` — scoring function for 7-point Likert
- **Ground truth:** ML-based classification from gibberish scale responses. The gibberish items are meaningless sentences rated on a Likert scale; respondents who give patterned responses to meaningless items are flagged.
- **ReReReRe suitability:** **Moderate.** Small N and few items in the gibberish data itself, but the Study1_IDs file may link to a larger survey dataset.

### 05_Bruhlmann_2020 — Data Quality in CrowdFlower
- **Paper:** Brühlmann, Petralito, Aeschbach & Opwis (2020). The quality of data collected online: An investigation of careless responding in a crowdsourced sample. *Methods in Psychology*.
- **Source:** https://osf.io/9vjur/
- **Files:**
  - `data_anon.csv` — 394 respondents, 163 columns (semicolon-delimited, CR line endings, Latin-1 encoding)
  - `Analysis.R` — analysis script
- **Ground truth:** Multiple indicators:
  - `quality` column — CrowdFlower platform quality rating
  - `v_IRI` — instructed response item (bogus item check)
  - `v_Bogus_Item` — another bogus item
  - `v_SRCR*`, `v_SRPR*`, `v_SRRR*`, `v_SRSI*` — self-reported carelessness scales
  - `v_Serious` — self-reported seriousness
  - `duration` — completion time
- **Scales:** PANAS (10 items), AttrakDiff (22 items), Need Fulfillment (24 items), Technology Commitment (14 items), Visual Aesthetics (14 items), BFI-44 (44 items), and more.
- **ReReReRe suitability:** **Excellent.** 163 items across many scales, multiple ground truth indicators, crowdsourced (high expected CR rate). This is probably the best Strategy 1 dataset.

### 06_Schneider_QoL — Careless Responding in Quality of Life Assessments
- **Paper:** Schneider, May & Stone. Careless responding in Internet-based quality of life assessments. *University of Southern California*.
- **Source:** https://osf.io/um9d3/
- **Files:**
  - `carersp.csv` — 1987 respondents, 58 columns
  - `variable_description.docx` — codebook
- **Scales:** PROMIS measures: Depression (8 items), Pain (10 items), Fatigue (7 items), Cognitive function (8 items). Plus demographics.
- **Ground truth:** Latent class analysis identified careless subgroups. Includes response time and vocabulary test scores as validity indicators.
- **ReReReRe suitability:** **Moderate.** Only 4 clinical domains (~33 items total), but interesting for clinical/health context. The latent class approach is a more sophisticated ground truth than simple attention checks.

### 07_Kuang_2025 — Individual Factors in Careless Responding (Mental Health)
- **Paper:** Kuang et al. (2025). The Impact of Individual Factors on Careless Responding Across Different Mental Disorder Screenings. *JMIR*.
- **Source:** https://zenodo.org/records/10423537
- **Files:** phq9.csv, pss.csv, isi.csv, gad7.csv, demographic.csv — each ~24,292 rows
- **Ground truth:** Page Time Index (PTI): items completed in < 2 seconds coded as careless.
- **ReReReRe suitability:** **Low.** Each scale is very short (7-14 items), and they are separate files (not a single multi-scale battery). Very large N though, and PTI provides item-level carelessness labels.

### 08_Niessen_2016 — Person-Fit Statistics for Careless Detection
- **Paper:** Niessen, Meijer & Tendeiro (2016). Detecting careless respondents in web-based questionnaires: Which method to use? *Journal of Research in Personality*.
- **Source:** https://dataverse.nl (doi:10.34894/F9ZCKM)
- **Files:** `Raw_data.sav` (SPSS), `Recoding_items.docx`
- **Ground truth:** Experimental manipulation — some respondents instructed to respond quickly.
- **ReReReRe suitability:** **Needs exploration.** SPSS format needs conversion. Paper is a classic in the field.

### gjkev (07_gjkev_CR_overtime) — Response Style Changes Over Time
- **Source:** https://osf.io/gjkev/
- **Files:** `Raw_data.sav`, `Clean_data.sav` (SPSS)
- **Design:** Studies changes in participants' responding styles over time.
- **ReReReRe suitability:** **Needs exploration.** Interesting longitudinal angle.

## Datasets NOT Downloaded (No Data Available or Simulations Only)

- **Curran (2016)** — OSF project exists (osf.io/6dkhm) but appears empty
- **Alfons & Welz (2024)** — GitHub repo is simulation-only (no empirical data)
- **Huang et al. (2012)** — Classic IER paper, data not publicly available
- **Meade & Craig (2012)** — Foundational paper, data not publicly shared
- **Biemann et al. (2025)** — Markov chain (Laz.R) method, data availability unclear
- **Kay (2024/2025)** — IDRIS/IDRIA validation, data may be on OSF but couldn't locate
- **Lazdauskas & McDevitt (2025)** — Temperament data, stated as available but couldn't access

## Recommended Priority for ReReReRe Validation

### Strategy 1 (Known ground truth):
1. **05_Bruhlmann_2020** — Best overall: 163 items, multiple ground truth types, crowdsourced
2. **01_Schroeders_2022** — Clean experimental labels, many items, widely cited
3. **03_Zoltak_2024** — Experimental + paradata, large sample
4. **08_Niessen_2016** — Classic paper, needs SPSS conversion

### Strategy 3 (Effect strengthening):
1. **06_Schneider_QoL** — Clinical scales with known psychometric properties
2. **05_Bruhlmann_2020** — BFI-44 factor structure as benchmark
3. **07_Kuang_2025** — Very large N, known clinical cutoffs

### Limitations to note:
- **Instructed carelessness ≠ natural carelessness:** Experimental induction (01, 03, 08) produces clean labels but may not capture realistic careless patterns
- **Small Ns in some:** Bloy (606), Brühlmann (394) limit statistical power
- **Few very-long instruments:** Most datasets have <200 items, which puts us below RR's sweet spot (nF≥15). Brühlmann (163 items) is closest.
