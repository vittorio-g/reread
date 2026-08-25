# Derive a homogeneous BOGUS/INSTRUCTED ground truth (>=5 bogus items) per dataset,
# aligned by row position to each dataset's scored matrix. Writes <name>_bogusgt.csv.
import csv, os
D = "../Dataset/gt_benchmark_candidates"
OUT = "bogus_bench"; os.makedirs(OUT, exist_ok=True)

def rdrows(path, **kw):
    with open(path, encoding="utf-8", errors="replace", newline="") as f:
        return list(csv.DictReader(f, **kw))

def num(v):
    try: return float(v)
    except: return None

def write_gt(name, gt, nmatrix, note):
    assert len(gt) == nmatrix, f"{name}: gt {len(gt)} != matrix {nmatrix}"
    with open(f"{OUT}/{name}_gt.csv", "w", newline="") as f:
        f.write("gt\n"); f.writelines(f"{g}\n" for g in gt)
    pos = sum(gt); print(f"{name:14s} n={len(gt):>6} pos={pos:>5} rate={pos/len(gt)*100:5.1f}%  [{note}]")

def matrix_n(path): return sum(1 for _ in open(path)) - 1

# ---- Kay S1/S2/S5: IDRIS 14-item infrequency, fail = endorsed (>=1); careless = >=2 fails ----
for k in (1, 2, 5):
    raw = rdrows(f"{D}/kay_idris_idria/study{k}_data.csv")
    icols = [h for h in raw[0].keys() if h.startswith("idris_") and h.replace("idris_", "").isdigit()]
    gt = []
    for r in raw:
        fails = sum(1 for c in icols if (num(r[c]) is not None and num(r[c]) >= 1))
        gt.append(1 if fails >= 2 else 0)
    nm = matrix_n(f"{D}/kay_idris_idria/kay_matrix_s{k}.csv")
    write_gt(f"kay_s{k}", gt, nm, f"IDRIS>=2 of {len(icols)}")

# ---- warning: 6 bogus (3 instructed iri + 3 infrequency if); careless = >=1 fail ----
raw = rdrows(f"{D}/warning_ipipneo300/data.csv")
correct = {"iri_1": 0, "iri_2": 2, "iri_3": 4}
gt = []
for r in raw:
    f = 0
    for c, corr in correct.items():
        v = num(r.get(c));
        if v is not None and v != corr: f += 1
    for c in ("if_1", "if_2", "if_3"):
        v = num(r.get(c))
        if v is not None and v >= 5: f += 1
    gt.append(1 if f >= 1 else 0)
write_gt("warning", gt, matrix_n(f"{D}/warning_ipipneo300/_matrix.csv"), "6 bogus (3 iri+3 if), fail>=1")

# ---- krause: existing 6-reactive-check label (fail>=2), built by feature_eng.py ----
lab = [int(float(r["y_ge2"])) for r in rdrows(f"{D}/krause_ier_youth/krause_labels.csv")]
write_gt("krause", lab, matrix_n(f"{D}/krause_ier_youth/krause_matrix.csv"), "6 reactive checks, fail>=2")

# ---- smarvus: existing 6-instructed label (fail>=1) ----
lab = [int(float(r["y"])) for r in rdrows(f"{D}/smarvus/smarvus_labels.csv")]
write_gt("smarvus", lab, matrix_n(f"{D}/smarvus/smarvus_matrix.csv"), "6 instructed, fail>=1")
print("done")
