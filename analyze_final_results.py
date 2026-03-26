"""Analyze Final Multiverse Results — ReReReRe + Benchmarks"""
import pandas as pd
import numpy as np
import os

os.chdir(os.path.dirname(os.path.abspath(__file__)))

# ============================================================
# 1. LOAD DATA
# ============================================================
print("Loading data...")
rr = pd.read_csv("final_multiverse_results.csv")
bench = pd.read_csv("final_benchmark_results.csv")
bytype = pd.read_csv("final_multiverse_results_by_type.csv")

print(f"ReReReRe results:  {len(rr):,} rows")
print(f"Benchmark results: {len(bench):,} rows")
print(f"By-type results:   {len(bytype):,} rows")
print(f"\nParameters used:")
for col in ['nFactors', 'n_respondents', 'pct_careless', 'corProp', 'z_threshold']:
    if col in rr.columns:
        vals = sorted(rr[col].unique())
        print(f"  {col}: {vals}")
print(f"  Replications: {rr['rep_id'].nunique()}")

# ============================================================
# 2. ReReReRe OVERALL SUMMARY
# ============================================================
print("\n" + "="*65)
print("  ReReReRe OVERALL SUMMARY")
print("="*65)

mcc = rr['mcc'].dropna()
print(f"  MCC:  mean={mcc.mean():.3f}  median={mcc.median():.3f}  max={mcc.max():.3f}")
print(f"  MCC >= 0.3: {(mcc >= 0.3).mean()*100:.1f}%")
print(f"  MCC >= 0.5: {(mcc >= 0.5).mean()*100:.1f}%")
print(f"  MCC >= 0.7: {(mcc >= 0.7).mean()*100:.1f}%")
print(f"  MCC <  0:   {(mcc < 0).mean()*100:.1f}%")

auc = rr.groupby(['rep_id','nFactors','n_respondents','pct_careless','corProp'])['auc'].first()
print(f"\n  AUC:  mean={auc.mean():.3f}  median={auc.median():.3f}")

# ============================================================
# 3. ReReReRe BY nFactors
# ============================================================
print("\n" + "-"*65)
print("  ReReReRe MCC by nFactors (averaged over all other params)")
print("-"*65)
by_nf = rr.groupby('nFactors')['mcc'].agg(['mean','median','std']).round(3)
print(by_nf.to_string())

# ============================================================
# 4. ReReReRe BY n_respondents
# ============================================================
print("\n" + "-"*65)
print("  ReReReRe MCC by n_respondents")
print("-"*65)
by_n = rr.groupby('n_respondents')['mcc'].agg(['mean','median','std']).round(3)
print(by_n.to_string())

# ============================================================
# 5. ReReReRe BY corProp
# ============================================================
print("\n" + "-"*65)
print("  ReReReRe MCC by corProp")
print("-"*65)
by_cp = rr.groupby('corProp')['mcc'].agg(['mean','median','std']).round(3)
print(by_cp.to_string())

# ============================================================
# 6. ReReReRe BY z_threshold
# ============================================================
print("\n" + "-"*65)
print("  ReReReRe MCC by z_threshold")
print("-"*65)
by_zt = rr.groupby('z_threshold')['mcc'].agg(['mean','median','std']).round(3)
print(by_zt.to_string())

# ============================================================
# 7. OPTIMAL z_threshold BY nFactors
# ============================================================
print("\n" + "-"*65)
print("  Optimal z_threshold by nFactors (highest mean MCC)")
print("-"*65)
nf_zt = rr.groupby(['nFactors','z_threshold'])['mcc'].mean().reset_index()
best_zt = nf_zt.loc[nf_zt.groupby('nFactors')['mcc'].idxmax()]
for _, row in best_zt.iterrows():
    print(f"  nF={int(row['nFactors']):2d} -> z_threshold={row['z_threshold']:.1f}  (mean MCC={row['mcc']:.3f})")

# ============================================================
# 8. SWEET SPOT (nF>=15, n>=300, corProp<=0.10, zt in [2.5,4])
# ============================================================
print("\n" + "-"*65)
print("  SWEET SPOT: nF>=15, n>=300, corProp<=0.10, z_threshold in [2.5, 4.0]")
print("-"*65)
sweet = rr[
    (rr['nFactors'] >= 15) &
    (rr['n_respondents'] >= 300) &
    (rr['corProp'] <= 0.10) &
    (rr['z_threshold'] >= 2.5) &
    (rr['z_threshold'] <= 4.0)
]
if len(sweet) > 0:
    print(f"  N cells: {len(sweet):,}")
    print(f"  MCC:         mean={sweet['mcc'].mean():.3f}  median={sweet['mcc'].median():.3f}")
    print(f"  Sensitivity: mean={sweet['sensitivity'].mean():.3f}")
    print(f"  Specificity: mean={sweet['specificity'].mean():.3f}")
else:
    print("  (no cells match — check parameter ranges)")

# ============================================================
# 9. RECOMMENDED ZONE (nF>=10, n>=100)
# ============================================================
print("\n" + "-"*65)
print("  RECOMMENDED ZONE: nF>=10, n>=100")
print("-"*65)
rec = rr[(rr['nFactors'] >= 10) & (rr['n_respondents'] >= 100)]
print(f"  N cells: {len(rec):,}")
print(f"  MCC: mean={rec['mcc'].mean():.3f}  median={rec['mcc'].median():.3f}")
print(f"  MCC >= 0.3: {(rec['mcc'] >= 0.3).mean()*100:.1f}%")
print(f"  MCC <  0:   {(rec['mcc'] < 0).mean()*100:.1f}%")

# ============================================================
# 10. BENCHMARK OVERALL SUMMARY
# ============================================================
print("\n" + "="*65)
print("  BENCHMARK COMPARISON (oracle best threshold)")
print("="*65)

bench_summary = bench.groupby('method').agg(
    mean_mcc=('best_mcc', 'mean'),
    median_mcc=('best_mcc', 'median'),
    mean_auc=('auc', 'mean'),
    median_auc=('auc', 'median'),
    mean_sens=('sensitivity', 'mean'),
    mean_spec=('specificity', 'mean')
).round(3)
print(bench_summary.to_string())

# ============================================================
# 11. BENCHMARK BY nFactors
# ============================================================
print("\n" + "-"*65)
print("  Benchmark MCC by nFactors")
print("-"*65)
bench_nf = bench.groupby(['nFactors','method'])['best_mcc'].mean().unstack('method').round(3)
print(bench_nf.to_string())

print("\n" + "-"*65)
print("  Benchmark AUC by nFactors")
print("-"*65)
bench_nf_auc = bench.groupby(['nFactors','method'])['auc'].mean().unstack('method').round(3)
print(bench_nf_auc.to_string())

# ============================================================
# 12. ReReReRe vs BENCHMARKS — HEAD-TO-HEAD BY nFactors
# ============================================================
print("\n" + "="*65)
print("  ReReReRe vs BENCHMARKS — HEAD-TO-HEAD")
print("="*65)

# Get ReReReRe best MCC per (rep, nF, n, pct) — optimize over corProp and z_threshold
rr_best = rr.groupby(['rep_id','nFactors','n_respondents','pct_careless'])['mcc'].max().reset_index()
rr_best.rename(columns={'mcc': 'best_mcc'}, inplace=True)
rr_best['method'] = 'ReReReRe'

# Also get ReReReRe AUC (best over corProp)
rr_auc_best = rr.groupby(['rep_id','nFactors','n_respondents','pct_careless','corProp'])['auc'].first().reset_index()
rr_auc_best = rr_auc_best.groupby(['rep_id','nFactors','n_respondents','pct_careless'])['auc'].max().reset_index()

rr_best = rr_best.merge(rr_auc_best, on=['rep_id','nFactors','n_respondents','pct_careless'])

# Combine with benchmark
bench_slim = bench[['rep_id','nFactors','n_respondents','pct_careless','method','best_mcc','auc']].copy()
combined = pd.concat([rr_best, bench_slim], ignore_index=True)

print("\nMean best MCC by nFactors x method:")
pivot_mcc = combined.groupby(['nFactors','method'])['best_mcc'].mean().unstack('method').round(3)
# Reorder columns
cols_order = ['ReReReRe','Mahalanobis','IRV','LongString','PersonTotal']
cols_order = [c for c in cols_order if c in pivot_mcc.columns]
pivot_mcc = pivot_mcc[cols_order]
print(pivot_mcc.to_string())

print("\nMean AUC by nFactors x method:")
pivot_auc = combined.groupby(['nFactors','method'])['auc'].mean().unstack('method').round(3)
pivot_auc = pivot_auc[[c for c in cols_order if c in pivot_auc.columns]]
print(pivot_auc.to_string())

# Win rates: ReReReRe vs Mahalanobis per cell
if 'Mahalanobis' in combined['method'].values:
    rr_cells = combined[combined['method']=='ReReReRe'][['rep_id','nFactors','n_respondents','pct_careless','best_mcc']].copy()
    mah_cells = combined[combined['method']=='Mahalanobis'][['rep_id','nFactors','n_respondents','pct_careless','best_mcc']].copy()
    merged = rr_cells.merge(mah_cells, on=['rep_id','nFactors','n_respondents','pct_careless'], suffixes=('_rr','_mah'))
    
    print("\n" + "-"*65)
    print("  ReReReRe vs Mahalanobis WIN RATE by nFactors")
    print("-"*65)
    for nf in sorted(merged['nFactors'].unique()):
        sub = merged[merged['nFactors']==nf]
        rr_wins = (sub['best_mcc_rr'] > sub['best_mcc_mah']).sum()
        mah_wins = (sub['best_mcc_mah'] > sub['best_mcc_rr']).sum()
        ties = len(sub) - rr_wins - mah_wins
        print(f"  nF={nf:2d}: ReReReRe wins {rr_wins}/{len(sub)} ({rr_wins/len(sub)*100:.0f}%), "
              f"Mah wins {mah_wins}/{len(sub)} ({mah_wins/len(sub)*100:.0f}%), ties {ties}")

# ============================================================
# 13. BY-TYPE DETECTION RATES
# ============================================================
print("\n" + "="*65)
print("  DETECTION BY CARELESS PATTERN TYPE")
print("="*65)

# Use a reasonable z_threshold range
zt_mid = bytype[(bytype['z_threshold'] >= 2.0) & (bytype['z_threshold'] <= 3.5)]
if len(zt_mid) > 0:
    print(f"\n  (z_threshold in [2.0, 3.5])")
    type_det = zt_mid.groupby('pattern')['detection_rate'].agg(['mean','median','std']).round(3)
    print(type_det.to_string())
    
    print("\n  Detection rate by pattern x nFactors:")
    # Need to merge nFactors from main results — bytype should have it
    if 'nFactors' in bytype.columns:
        type_nf = zt_mid.groupby(['nFactors','pattern'])['detection_rate'].mean().unstack('pattern').round(3)
        print(type_nf.to_string())

# ============================================================
# 14. REPLICATION STABILITY
# ============================================================
print("\n" + "="*65)
print("  REPLICATION STABILITY")
print("="*65)

# SD of MCC across replications for same condition
rep_sd = rr.groupby(['nFactors','n_respondents','pct_careless','corProp','z_threshold'])['mcc'].std().dropna()
print(f"  Mean within-condition SD of MCC: {rep_sd.mean():.3f}")
print(f"  Median within-condition SD:      {rep_sd.median():.3f}")
print(f"  95th percentile SD:              {rep_sd.quantile(0.95):.3f}")
print(f"  SE of mean (R=50):               ~{rep_sd.mean()/np.sqrt(50):.4f}")

# ============================================================
# 15. RECOMMENDED ZONE WITH BEST corProp
# ============================================================
print("\n" + "-"*65)
print("  Best corProp overall and by nFactors")
print("-"*65)
cp_mcc = rr.groupby('corProp')['mcc'].mean()
print(f"  Overall best corProp: {cp_mcc.idxmax()} (mean MCC={cp_mcc.max():.3f})")

cp_by_nf = rr.groupby(['nFactors','corProp'])['mcc'].mean().reset_index()
best_cp = cp_by_nf.loc[cp_by_nf.groupby('nFactors')['mcc'].idxmax()]
for _, row in best_cp.iterrows():
    print(f"  nF={int(row['nFactors']):2d} -> best corProp={row['corProp']:.2f}  (mean MCC={row['mcc']:.3f})")

print("\n" + "="*65)
print("  ANALYSIS COMPLETE")
print("="*65)
