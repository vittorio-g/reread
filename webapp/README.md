# CaReReRe web tool — re-re.re

**CaReReRe** = *Careless Respondent Recognition via Resampled Ensemble* (brand name;
the core permutation index is `rr`). Client-side web app that flags careless
respondents in questionnaire data using the ensemble (`rr` + longstring +
person-total, logistic combiner). **All computation runs in the visitor's browser
— no data is ever transmitted.** (JS namespace stays `ReReRe` internally.)

## Layout

```
webapp/
├── site/               <- the deployable static site (4 files, no build step)
│   ├── index.html      UI + styles
│   ├── app.js          UI logic (drag&drop, table, histogram, downloads)
│   ├── worker.js       Web Worker wrapper (keeps UI responsive)
│   └── rerere.js       core algorithm library (browser + Node)
├── test_core.js        Node test-suite (synthetic + fidelity vs R pipeline)
└── README.md
```

## Validation

`node test_core.js` — 10/10:
- AUC 0.955 on synthetic factor data with 20% injected random careless
- noise-calibrated structure diagnostic (warns on structureless data)
- CSV parser: quoted fields, `,`/`;`/tab sniffing, EU decimals
- **fidelity vs the R pipeline on the real study data (N=84, 109 items):
  Spearman rho = 0.984**

## Deploy (Cloudflare Pages)

One-time setup (user):
```bash
npm i -g wrangler
wrangler login          # browser OAuth on the Cloudflare account holding re-re.re
```

Deploy (repeatable):
```bash
cd webapp
wrangler pages project create rerere --production-branch=main   # first time only
wrangler pages deploy site --project-name=rerere
```

Custom domain (first time only): Cloudflare dashboard → Workers & Pages →
`rerere` → Custom domains → add `re-re.re` (one click — the domain is already
on Cloudflare, DNS is handled automatically).

## Algorithm defaults (paper-canonical)

corProp 0.03 · permutation iterations 200 · min pairs 15 · z threshold 1.5
(low z = careless) · proportion rescale · automatic sign-alignment of
reverse-coded items · median imputation for missing (n_missing reported).

Advanced panel exposes: corProp, threshold (live re-flag without recompute),
iterations, min pairs, seed, delimiter override.

## Roadmap

- XLSX input (SheetJS) — currently "save as CSV"
- ensemble view (z_RR + auxiliaries with per-index thresholds)
- R package (`ReReReRe` on CRAN)
