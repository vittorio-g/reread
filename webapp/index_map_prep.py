"""Prep real scale keys + aligned response-time for the index map.
Writes webapp/keys/<name>_scale.csv (one scale-id per matrix column, -1 = excluded)
and webapp/keys/<name>_rt.csv (one RT in seconds per matrix row, full length).
RT alignment is reproduced by the complete-case mask on the item columns and VERIFIED
by comparing matrix values to the raw-masked values before writing."""
import pandas as pd, numpy as np, json, os, re
os.makedirs("keys", exist_ok=True)
D = "../Dataset/gt_benchmark_candidates"

def matrix_header(path):
    return pd.read_csv(path, nrows=0).columns.tolist()

def write_scale(name, scale_ids):
    pd.DataFrame({"scale": scale_ids}).to_csv(f"keys/{name}_scale.csv", index=False)
    k = len({s for s in scale_ids if s >= 0})
    print(f"  scale {name}: {len(scale_ids)} items -> {k} scales, {sum(1 for s in scale_ids if s<0)} excluded")

def verify_and_write_rt(name, mat_path, raw, item_cols_raw, mat_cols, rt_col):
    """mask raw to complete-case on item_cols_raw (0->NaN), reorder to mat_cols, verify == matrix."""
    m = raw[item_cols_raw].replace(0, np.nan)
    mask = m.dropna().index
    sub = raw.loc[mask]
    M = pd.read_csv(mat_path)
    if len(sub) != len(M):
        print(f"  ! RT {name}: masked {len(sub)} != matrix {len(M)} — SKIP"); return
    # map matrix col -> raw col (same names)
    common = [c for c in mat_cols if c in raw.columns]
    a = sub[common].reset_index(drop=True).astype(float).values
    b = M[common].reset_index(drop=True).astype(float).values
    if not np.allclose(a, b, equal_nan=True, atol=1e-6):
        diff = np.nanmean(np.abs(a - b)); print(f"  ! RT {name}: values differ (mean|Δ|={diff:.3f}) — SKIP"); return
    rt = pd.to_numeric(sub[rt_col].reset_index(drop=True), errors="coerce")
    pd.DataFrame({"rt": rt}).to_csv(f"keys/{name}_rt.csv", index=False)
    print(f"  RT {name}: aligned {len(rt)} rows (verified), median={np.nanmedian(rt):.0f}s")

# ---------- Study1: real 11 domains + reverse from items.json ----------
it = json.load(open("../Raccolta Dati/Careless_Responding_Study/docs/items.json", encoding="utf-8"))["items"]
dom = {x["id"]: x["domain"] for x in it}
cols = [c.strip('"') for c in matrix_header("_study_matrix.csv")]
doms = sorted({dom[c] for c in cols if dom.get(c) and dom[c] != "Random"})
dmap = {d: i for i, d in enumerate(doms)}
write_scale("Study1_induced", [dmap.get(dom.get(c, "Random"), -1) if dom.get(c) != "Random" else -1 for c in cols])

# ---------- 16pf: letter prefix = factor ----------
cols = matrix_header(f"{D}/opsy_16pf_matrix.csv")
letters = sorted({re.match(r"^([A-P])\d+$", c).group(1) for c in cols if re.match(r"^([A-P])\d+$", c)})
lmap = {l: i for i, l in enumerate(letters)}
write_scale("opsy_16pf", [lmap.get(re.match(r"^([A-P])\d+$", c).group(1), -1) if re.match(r"^([A-P])\d+$", c) else -1 for c in cols])

# ---------- duckworth: E/N/A/C/O prefix ----------
cols = matrix_header(f"{D}/duckworth_grit_vcl/duckworth_matrix.csv")
tr = sorted({c[0] for c in cols if c[0] in "ENACO"})
tmap = {t: i for i, t in enumerate(tr)}
write_scale("duckworth", [tmap.get(c[0], -1) if c[0] in "ENACO" else -1 for c in cols])

# ---------- hexaco: facet prefix (letters before trailing digits) ----------
cols = matrix_header(f"{D}/opsy_hexaco_matrix.csv")
fac = lambda c: re.match(r"^(.*?)\d+$", c).group(1) if re.match(r"^(.*?)\d+$", c) else None
facets = sorted({fac(c) for c in cols if fac(c)})
fmap = {f: i for i, f in enumerate(facets)}
write_scale("opsy_hexaco", [fmap.get(fac(c), -1) for c in cols])

print("\n-- RT (verified alignment) --")
# duckworth: testelapse
raw = pd.read_csv(f"{D}/duckworth_grit_vcl/data.csv", sep="\t", low_memory=False)
ip = [c for c in matrix_header(f"{D}/duckworth_grit_vcl/duckworth_matrix.csv")]
verify_and_write_rt("duckworth", f"{D}/duckworth_grit_vcl/duckworth_matrix.csv", raw, ip, ip, "testelapse")
# 16pf: elapsed
raw = pd.read_csv(f"{D}/opsy_16pf/data.csv", sep="\t", low_memory=False)
ip = matrix_header(f"{D}/opsy_16pf_matrix.csv")
verify_and_write_rt("opsy_16pf", f"{D}/opsy_16pf_matrix.csv", raw, ip, ip, "elapsed")
# hexaco: elapse
raw = pd.read_csv(f"{D}/opsy_hexaco/data.csv", sep="\t", low_memory=False)
ip = matrix_header(f"{D}/opsy_hexaco_matrix.csv")
verify_and_write_rt("opsy_hexaco", f"{D}/opsy_hexaco_matrix.csv", raw, ip, ip, "elapse")
# warning: duration
raw = pd.read_csv(f"{D}/warning_ipipneo300/data.csv", low_memory=False)
ip = matrix_header(f"{D}/warning_ipipneo300/_matrix.csv")
verify_and_write_rt("warning_IPIP300", f"{D}/warning_ipipneo300/_matrix.csv", raw, ip, ip, "duration")
print("done")
