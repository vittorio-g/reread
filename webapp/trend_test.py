# trend_test.py — ONE p-value for the strictness thesis (report roadmap item 6).
# Pooled row-level model on the healthy datasets (smarvus, Kay S1), estimable rows, band <=25%:
#   piHat = a_d + b_d * true + gamma * (k - kbar_d) * true
# gamma = common change of the injection slope per unit of criterion strictness.
# H0: gamma <= 0, tested by cluster bootstrap over targets within (dataset, k) — the honest
# resampling unit (resamples within a target share pools).
import csv, random, collections
random.seed(7)
rows = [r for r in csv.DictReader(open("strictness_sweep.csv"))
        if r["dataset"] in ("smarvus", "Kay S1") and r["estimable"] == "1"
        and float(r["true_prev"]) <= 0.25]
for r in rows:
    r["x"] = float(r["true_prev"]); r["y"] = float(r["piHat"]); r["k"] = int(r["k"])
kbar = {d: sum(r["k"] for r in rows if r["dataset"] == d) / sum(1 for r in rows if r["dataset"] == d)
        for d in ("smarvus", "Kay S1")}

def fit(rs):
    # design: [1_sm, 1_ks, x_sm, x_ks, x*(k-kbar_d)] — solve normal equations (5x5)
    X, Y = [], []
    for r in rs:
        sm = 1.0 if r["dataset"] == "smarvus" else 0.0
        X.append([sm, 1 - sm, r["x"] * sm, r["x"] * (1 - sm),
                  r["x"] * (r["k"] - kbar[r["dataset"]])])
        Y.append(r["y"])
    p = 5
    A = [[sum(X[i][a] * X[i][b] for i in range(len(X))) for b in range(p)] for a in range(p)]
    B = [sum(X[i][a] * Y[i] for i in range(len(X))) for a in range(p)]
    for c in range(p):                      # gaussian elimination with partial pivot
        piv = max(range(c, p), key=lambda rr: abs(A[rr][c]))
        A[c], A[piv] = A[piv], A[c]; B[c], B[piv] = B[piv], B[c]
        if abs(A[c][c]) < 1e-12: return None
        for rr in range(c + 1, p):
            f = A[rr][c] / A[c][c]
            for cc in range(c, p): A[rr][cc] -= f * A[c][cc]
            B[rr] -= f * B[c]
    beta = [0.0] * p
    for c in range(p - 1, -1, -1):
        beta[c] = (B[c] - sum(A[c][cc] * beta[cc] for cc in range(c + 1, p))) / A[c][c]
    return beta

beta = fit(rows)
print("pooled model rows:", len(rows))
print("base slopes: smarvus %.3f, Kay S1 %.3f" % (beta[2], beta[3]))
print("gamma (d slope / d k): %+.4f" % beta[4])

clusters = collections.defaultdict(list)
for r in rows: clusters[(r["dataset"], r["k"], r["x"])].append(r)
by_dk = collections.defaultdict(list)
for (d, k, x) in clusters: by_dk[(d, k)].append((d, k, x))
gs = []
for _ in range(4000):
    bs = []
    for dk, keys in by_dk.items():
        for _ in keys: bs.extend(clusters[keys[random.randrange(len(keys))]])
    b = fit(bs)
    if b: gs.append(b[4])
gs.sort()
lo, hi = gs[int(.025 * len(gs))], gs[int(.975 * len(gs))]
p_one = sum(1 for g in gs if g <= 0) / len(gs)
print("cluster-bootstrap 95%% CI: [%+.4f, %+.4f]   one-sided p (gamma<=0): %.4f" % (lo, hi, p_one))
