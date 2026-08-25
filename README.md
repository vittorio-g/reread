# reread — screening careless respondents with a self-calibrating threshold

Research repository for *Beyond Ranking Careless Respondents: A Self-Calibrating
Threshold and a Structure-Free Index for Screening Questionnaire Data*
(Guerrieri, Gallucci & Passarelli).

Careless screening asks the analyst three questions that have no default answer:
which index to compute, how to combine indices, and where to cut. This project
closes all three. The threshold is anchored on the attentive mode of the score
distribution — estimated from that mode's left flank alone, where the careless
cannot reach — and is reported together with a calibrated prevalence estimate and
with what the same rule flags on careless-free data of the same size. The score it
cuts is `rc` (*relative coherence*), a permutation-based index of within-person
coherence that needs no item-to-factor map, combined with LongString and
Person–Total in a fixed-weight, profile-gated ensemble.

## Layout

| Path | What it is |
|---|---|
| `reread/` | The R package. `reread()` is the one-call entry point; `auto_flag()` is the label-free calibration; `rc_index()` the index alone. |
| `webapp/site/` | The browser tool served at <https://ca.re-re.re> — `rerere.js` is the reference engine, and everything runs client-side. |
| `webapp/` | Analysis scripts for the three studies, the external benchmarks, and the simulation. |
| `thesis/`, `paper_assets_v3/`, `article_assets_v3/` | Figures and assets for the manuscript and thesis. |
| `archive/` | Superseded simulation runs and reports, kept for the record. |
| `CLAUDE.md` | The working log: every experiment run, including the ones that failed and why. |

## Reproducing

```r
# from a clone
install.packages("reread_0.1.0.tar.gz", repos = NULL, type = "source")
library(reread)
fit <- reread(your_responses)      # scores, flags, calibrated prevalence
fit$scores
```

The R package and the browser engine agree to six decimal places on the attentive
mode, the scale, the cut, the prevalence estimate and the flag set; the parity
check lives in `webapp/`.

## Data not in this repository

Third-party benchmark corpora (`Dataset/`) and the item-level matrices derived
from them are excluded: they are large, redistributable only under their own
licences, and re-downloadable from the sources named in each subdirectory README.
The derived per-respondent scores that the analyses actually read are committed
under `webapp/ext_bench2/*_y.csv`.

This repository is **private**. It contains item-level responses from the
validation study, which is embargoed until the manuscript is accepted; the public
release is prepared separately on OSF.
