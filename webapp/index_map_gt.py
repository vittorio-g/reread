"""Interpret the index-space factors against INDEPENDENT ground truth.
Per respondent: factor scores = standardized 9 indices @ varimax loadings.
Spearman-correlate each factor with each dataset's external criterion, grouped by
criterion type (induced / infrequency-bogus / attention-check / overclaiming / speeding-RT).
Factor scores oriented high=careless; positive r => that factor tracks that criterion."""
import numpy as np, pandas as pd
from scipy.stats import spearmanr

COLS = ["rc","synonyms","antonyms","even_odd","rpr","person_total","d2","longstring","irv"]
long = pd.read_csv("index_map_long.csv")
long[COLS + ["response_time"]] = long[COLS + ["response_time"]].apply(pd.to_numeric, errors="coerce")
Ld = pd.read_csv("index_map_loadings.csv", index_col=0)
FAC = list(Ld.columns)  # F1:pair-consistency, F2:normativity, F3:outlier/variance, F4:split-half
L = Ld.loc[COLS].values

D = "../Dataset/gt_benchmark_candidates"
# dataset -> list of (criterion_type, path, column, orient)  orient: +1 careless=high, -1 careless=low(e.g. time)
GT = {
  "Study1_induced":  [("induced",         "_study_labels.csv", "careless", +1)],
  "Kay_S1":          [("infrequency",     "bogus_bench/kay_s1_gt.csv", "gt", +1)],
  "Kay_S2":          [("infrequency",     "bogus_bench/kay_s2_gt.csv", "gt", +1)],
  "Kay_S5":          [("infrequency",     "bogus_bench/kay_s5_gt.csv", "gt", +1)],
  "warning_IPIP300": [("infrequency",     "bogus_bench/warning_gt.csv", "gt", +1)],
  "krause":          [("attention_check", "bogus_bench/krause_gt.csv", "gt", +1)],
  "douglas":         [("attention_check", D+"/douglas2023_dataquality/douglas_labels.csv", "y_attn_ge2", +1),
                      ("speeding_RT",     D+"/douglas2023_dataquality/douglas_labels.csv", "y_speeder", +1)],
  "duckworth":       [("overclaiming",    D+"/duckworth_grit_vcl/duckworth_gt.csv", "fake_ge2", +1)],
}
# continuous RT criterion straight from the aligned response_time column (oriented high=careless=fast)
RT_DATASETS = ["duckworth", "opsy_16pf", "warning_IPIP300"]

def fscore(sub):
    Z = sub[COLS].copy()
    Z = (Z - Z.mean()) / Z.std(ddof=0)     # standardize each index within dataset
    Z = Z.fillna(0.0)                       # mean-impute missing (e.g. antonyms NaN)
    return Z.values @ L                     # n x nFactors

rows = []
print("dataset          criterion         n   rate   " + "  ".join(f.split(':')[1][:9].ljust(9) for f in FAC))
for ds, crits in GT.items():
    sub = long[long.dataset == ds].reset_index(drop=True)
    if len(sub) == 0: continue
    F = fscore(sub)
    for ctype, path, col, orient in crits:
        try:
            g = pd.to_numeric(pd.read_csv(path)[col], errors="coerce").values
        except Exception as e:
            print(f"  ! {ds}/{col}: {e}"); continue
        if len(g) != len(sub):
            print(f"  ! {ds}/{col}: len {len(g)} != {len(sub)} (skip)"); continue
        gg = orient * g
        rr = []
        for f in range(len(FAC)):
            ok = np.isfinite(F[:, f]) & np.isfinite(gg)
            rho = spearmanr(F[ok, f], gg[ok]).correlation if ok.sum() > 10 else np.nan
            rr.append(rho)
        rate = np.nanmean((g > np.nanmedian(g)) if orient < 0 else (g > 0))
        rows.append([ds, ctype, len(sub), rate] + rr)
        print(f"{ds:<16} {ctype:<15} {len(sub):>5} {rate:5.2f}   " + "  ".join(f"{x:+.2f}".ljust(9) if np.isfinite(x) else "  NA     " for x in rr))
    # continuous RT criterion from the aligned response_time column
    if ds in RT_DATASETS and sub["response_time"].notna().sum() > 20:
        gg = sub["response_time"].values  # already oriented high=careless=fast
        rr = []
        for f in range(len(FAC)):
            ok = np.isfinite(F[:, f]) & np.isfinite(gg)
            rr.append(spearmanr(F[ok, f], gg[ok]).correlation if ok.sum() > 10 else np.nan)
        rows.append([ds, "speeding_RT", len(sub), np.nan] + rr)
        print(f"{ds:<16} {'speeding_RT':<15} {len(sub):>5}   -    " + "  ".join(f"{x:+.2f}".ljust(9) if np.isfinite(x) else "  NA     " for x in rr))

R = pd.DataFrame(rows, columns=["dataset","criterion","n","rate"]+FAC)
R.to_csv("index_map_gt.csv", index=False)

# pooled by criterion type (Fisher-z, unweighted)
print("\n=== Pooled by criterion type (Fisher-z mean) ===")
print("criterion         k  " + "  ".join(f.split(':')[1][:9].ljust(9) for f in FAC))
fisher=lambda r: np.arctanh(np.clip(r,-0.999,0.999))
for ct in ["induced","infrequency","attention_check","overclaiming","speeding_RT"]:
    d = R[R.criterion == ct]
    if len(d) == 0: continue
    means = [np.tanh(np.nanmean(fisher(d[f].values))) for f in FAC]
    print(f"{ct:<15} {len(d):>2}  " + "  ".join(f"{x:+.2f}".ljust(9) for x in means))
print("\nfactors:", ", ".join(FAC))
