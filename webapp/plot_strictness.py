# plot_strictness.py — analysis + figures for the bogus-strictness study.
# Metric: OLS slope of pi-hat on true injected prevalence, per (dataset, k), computed on
# ESTIMABLE resamples (declared-failure resamples are counted separately as availability).
# CI: 95% percentile bootstrap over resample rows (2000 draws). Companion: Pearson r.
# Fig A: slope +/- CI versus k, one series per dataset (slope 1 = ideal responsiveness).
# Fig B: raw injection clouds for smarvus, one facet per k, with the fitted line.
import csv, random, collections
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
plt.rcParams.update({"font.family": "serif", "font.size": 10.5})
random.seed(42)

rows = list(csv.DictReader(open("strictness_sweep.csv")))
for r in rows:
    r["true"] = float(r["true_prev"]); r["pi"] = float(r["piHat"])
    r["est"] = r["estimable"] == "1"; r["k"] = int(r["k"])

def ols_slope(pts):
    n = len(pts)
    mx = sum(p[0] for p in pts) / n; my = sum(p[1] for p in pts) / n
    sxx = sum((p[0]-mx)**2 for p in pts)
    if sxx == 0: return None
    return sum((p[0]-mx)*(p[1]-my) for p in pts) / sxx

def pearson(pts):
    n = len(pts); mx = sum(p[0] for p in pts)/n; my = sum(p[1] for p in pts)/n
    sxx = sum((p[0]-mx)**2 for p in pts); syy = sum((p[1]-my)**2 for p in pts)
    if sxx == 0 or syy == 0: return float("nan")
    return sum((p[0]-mx)*(p[1]-my) for p in pts) / (sxx*syy) ** 0.5

groups = collections.OrderedDict()
for r in rows:
    groups.setdefault((r["dataset"], r["k"]), []).append(r)

res = []
print("%-14s k  rows  estim%%   slope   [95%% CI]        r" % "dataset")
for (ds, k), g in groups.items():
    est = [(r["true"], r["pi"]) for r in g if r["est"]]
    frac = len(est) / len(g)
    if len(est) < 8 or len(set(p[0] for p in est)) < 3:
        print("%-14s %d  %4d  %5.0f%%   (too few estimable rows)" % (ds, k, len(g), 100*frac))
        res.append((ds, k, frac, None, None, None, None)); continue
    sl = ols_slope(est); rr = pearson(est)
    boots = []
    for _ in range(2000):
        bs = [est[random.randrange(len(est))] for _ in range(len(est))]
        b = ols_slope(bs)
        if b is not None: boots.append(b)
    boots.sort()
    lo, hi = boots[int(.025*len(boots))], boots[int(.975*len(boots))]
    res.append((ds, k, frac, sl, lo, hi, rr))
    print("%-14s %d  %4d  %5.0f%%   %+.3f  [%+.3f, %+.3f]  %.3f" % (ds, k, len(g), 100*frac, sl, lo, hi, rr))

with open("strictness_slopes.csv", "w", newline="") as f:
    w = csv.writer(f); w.writerow(["dataset","k","estimable_frac","slope","ci_lo","ci_hi","pearson_r"])
    for row in res: w.writerow(row)

# ---------------- Fig A: slope vs k ----------------
SERIES = [("warning IPIP", "#B35806", "o"), ("smarvus", "#2166AC", "s"),
          ("Kay S1", "#762A83", "^"), ("Kay S2", "#1B7837", "D")]
fig, ax = plt.subplots(figsize=(8.6, 5.2))
ax.axhline(1.0, color="#8a8a8a", lw=1.2, ls="--"); ax.text(0.62, 1.02, "ideal (slope 1)", color="#8a8a8a", fontsize=8.5)
ax.axhline(0.0, color="#cccccc", lw=1.0)
for ds, col, mk in SERIES:
    pts = [(k, sl, lo, hi) for (d, k, fr, sl, lo, hi, r) in res if d == ds and sl is not None]
    if not pts: continue
    xs = [p[0] for p in pts]; ys = [p[1] for p in pts]
    ax.errorbar(xs, ys, yerr=[[y-l for _, y, l, _ in pts], [h-y for _, y, _, h in pts]],
                color=col, marker=mk, ms=7, lw=2, capsize=4, label=ds, zorder=3)
    ax.annotate(ds, xy=(xs[-1], ys[-1]), xytext=(8, 0), textcoords="offset points",
                color=col, fontsize=9, va="center", fontweight="bold")
ax.set_xlabel("Carelessness criterion: minimum number of failed checks ($k$)")
ax.set_ylabel(r"Slope of $\hat\pi$ on injected careless fraction")
ax.set_title("Stricter ground truth, sharper estimate:\nresponsiveness of the prevalence estimate by criterion severity", fontsize=11.5)
ax.set_xticks([1, 2, 3, 4, 5]); ax.grid(alpha=.25); ax.legend(fontsize=9, loc="lower right")
plt.tight_layout(); plt.savefig("strictness_slopes.png", dpi=160); plt.close()

# ---------------- Fig B: smarvus clouds by k ----------------
ks = sorted({k for (d, k) in groups if d == "smarvus"})
fig, axes = plt.subplots(1, len(ks), figsize=(2.9*len(ks), 3.4), sharex=True, sharey=True)
for ax, k in zip(axes, ks):
    g = groups[("smarvus", k)]
    est = [(r["true"], r["pi"]) for r in g if r["est"]]
    non = [(r["true"], r["pi"]) for r in g if not r["est"]]
    ax.plot([0, .45], [0, .45], ls="--", color="#aaaaaa", lw=1)
    if non: ax.scatter([p[0] for p in non], [p[1] for p in non], s=14, facecolors="none",
                       edgecolors="#999999", linewidths=.8, label="non estimable")
    ax.scatter([p[0] for p in est], [p[1] for p in est], s=14, color="#2166AC", alpha=.65, label="estimable")
    sl = ols_slope(est) if len(est) >= 8 else None
    if sl is not None:
        mx = sum(p[0] for p in est)/len(est); my = sum(p[1] for p in est)/len(est)
        xs = [0, .45]; ax.plot(xs, [my + sl*(x-mx) for x in xs], color="#B35806", lw=2)
        ax.set_title("$k\\geq%d$   slope %.2f" % (k, sl), fontsize=10)
    else:
        ax.set_title("$k\\geq%d$" % k, fontsize=10)
    ax.grid(alpha=.2); ax.set_xlim(0, .45); ax.set_ylim(0, .55)
    ax.set_xlabel("true injected")
axes[0].set_ylabel(r"$\hat\pi$")
axes[0].legend(fontsize=7.5, loc="upper left")
fig.suptitle("smarvus: injection clouds by criterion severity (every resample shown)", fontsize=11)
plt.tight_layout(); plt.savefig("strictness_clouds.png", dpi=160, bbox_inches="tight"); plt.close()
print("wrote strictness_slopes.png, strictness_clouds.png, strictness_slopes.csv")
