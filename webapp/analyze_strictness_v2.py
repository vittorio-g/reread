# analyze_strictness_v2.py — full quality dashboard for the strictness study, replacing the
# plain-OLS analysis. Per (dataset, k):
#   slope_band   : WLS slope (weights = sample size n) on ESTIMABLE resamples with true <= 25%
#                  (the realistic band; keeps the heavy-contamination collapse regime out)
#   ci           : 95% CLUSTER bootstrap by target (resamples within a target share pools, so
#                  row bootstrap understates uncertainty; clusters are the honest unit)
#   sigma_resid  : weighted SD of residuals after the affine fit, in prevalence points
#                  (precision of a single estimate once the dataset's own baseline is absorbed)
#   theil_sen    : robust slope (median of pairwise slopes) as a robustness column
#   slope_global : WLS slope over the full 2-40% range (sensitivity)
#   avail        : share of resamples with an estimable pi-hat (reported, never mixed in)
import csv, random, collections
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
plt.rcParams.update({"font.family": "serif", "font.size": 10.5})
random.seed(42)
BAND = 0.25

rows = list(csv.DictReader(open("strictness_sweep.csv")))
for r in rows:
    r["true"] = float(r["true_prev"]); r["pi"] = float(r["piHat"])
    r["est"] = r["estimable"] == "1"; r["k"] = int(r["k"]); r["n"] = int(r["n"])

def wls(pts):  # pts = (x, y, w) -> (slope, intercept) or None
    sw = sum(p[2] for p in pts)
    mx = sum(p[0]*p[2] for p in pts)/sw; my = sum(p[1]*p[2] for p in pts)/sw
    sxx = sum(p[2]*(p[0]-mx)**2 for p in pts)
    if sxx <= 0: return None
    b = sum(p[2]*(p[0]-mx)*(p[1]-my) for p in pts)/sxx
    return b, my - b*mx

def resid_sd(pts, fit):
    b, a = fit; sw = sum(p[2] for p in pts)
    return (sum(p[2]*(p[1]-(a+b*p[0]))**2 for p in pts)/sw) ** 0.5

def theil_sen(pts):
    sl = []
    for i in range(len(pts)):
        for j in range(i+1, len(pts)):
            dx = pts[j][0]-pts[i][0]
            if abs(dx) > 1e-9: sl.append((pts[j][1]-pts[i][1])/dx)
    if not sl: return None
    sl.sort(); return sl[len(sl)//2]

groups = collections.OrderedDict()
for r in rows: groups.setdefault((r["dataset"], r["k"]), []).append(r)

res = []
hdr = "%-13s k  avail%% | slope(0-25) [cluster 95%% CI]  sigResid  TheilSen | slope(all)"
print(hdr % "dataset")
for (ds, k), g in groups.items():
    avail = sum(1 for r in g if r["est"]) / len(g)
    band = [(r["true"], r["pi"], r["n"]) for r in g if r["est"] and r["true"] <= BAND]
    allr = [(r["true"], r["pi"], r["n"]) for r in g if r["est"]]
    tb = sorted(set(p[0] for p in band))
    if len(band) < 8 or len(tb) < 3:
        print("%-13s %d  %5.0f%% | (banda insufficiente: %d righe, %d target)" % (ds, k, 100*avail, len(band), len(tb)))
        res.append(dict(ds=ds, k=k, avail=avail, slope=None)); continue
    fit = wls(band); slope = fit[0]
    sres = resid_sd(band, fit)
    ts = theil_sen(band)
    gfit = wls(allr); gslope = gfit[0] if gfit else float("nan")
    # cluster bootstrap by target
    by_t = collections.defaultdict(list)
    for p in band: by_t[p[0]].append(p)
    keys = list(by_t); boots = []
    for _ in range(2000):
        sel = [keys[random.randrange(len(keys))] for _ in keys]
        pts = [p for t in sel for p in by_t[t]]
        if len(set(p[0] for p in pts)) < 2: continue
        f = wls(pts)
        if f: boots.append(f[0])
    boots.sort()
    lo, hi = boots[int(.025*len(boots))], boots[int(.975*len(boots))]
    print("%-13s %d  %5.0f%% |   %+.3f    [%+.3f, %+.3f]    %.3f     %+.3f  |   %+.3f" %
          (ds, k, 100*avail, slope, lo, hi, sres, ts, gslope))
    res.append(dict(ds=ds, k=k, avail=avail, slope=slope, lo=lo, hi=hi,
                    sres=sres, ts=ts, gslope=gslope))

with open("strictness_dashboard.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["dataset","k","availability","slope_band","ci_lo","ci_hi","sigma_resid","theil_sen","slope_global"])
    for r in res:
        w.writerow([r["ds"], r["k"], "%.3f" % r["avail"]] +
                   (["","","","","",""] if r["slope"] is None else
                    ["%.4f" % r["slope"], "%.4f" % r["lo"], "%.4f" % r["hi"],
                     "%.4f" % r["sres"], "%.4f" % r["ts"], "%.4f" % r["gslope"]]))

# ---------------- dashboard figure: 3 stacked panels ----------------
SERIES = [("warning IPIP", "#B35806", "o"), ("smarvus", "#2166AC", "s"),
          ("Kay S1", "#762A83", "^"), ("Kay S2", "#1B7837", "D")]
fig, (axS, axR, axA) = plt.subplots(3, 1, figsize=(8.8, 9.2), sharex=True,
                                    gridspec_kw={"height_ratios": [3, 1.6, 1.2]})
axS.axhline(1.0, color="#8a8a8a", lw=1.2, ls="--")
axS.text(4.55, 1.03, "ideal (slope 1)", color="#8a8a8a", fontsize=8.5, ha="right")
axS.axhline(0.0, color="#cccccc", lw=1.0)
for ds, col, mk in SERIES:
    pts = [r for r in res if r["ds"] == ds and r["slope"] is not None]
    if not pts: continue
    xs = [r["k"] for r in pts]
    axS.errorbar(xs, [r["slope"] for r in pts],
                 yerr=[[r["slope"]-r["lo"] for r in pts], [r["hi"]-r["slope"] for r in pts]],
                 color=col, marker=mk, ms=7, lw=2, capsize=4, label=ds, zorder=3)
    axS.plot(xs, [r["ts"] for r in pts], color=col, marker=mk, ms=4, lw=0, alpha=.45,
             markerfacecolor="none", zorder=2)   # Theil-Sen as hollow companion
    axR.plot(xs, [100*r["sres"] for r in pts], color=col, marker=mk, ms=6, lw=1.8)
    av = [r for r in res if r["ds"] == ds]
    axA.plot([r["k"] for r in av], [100*r["avail"] for r in av], color=col, marker=mk, ms=6, lw=1.8)
axS.set_ylabel("Slope of $\\hat\\pi$ (band 0–25%)\nfilled = WLS, hollow = Theil–Sen")
axS.set_title("Prevalence-estimate quality dashboard by criterion severity\n"
              "(WLS, cluster bootstrap by target)", fontsize=11.5)
axS.grid(alpha=.25); axS.legend(fontsize=9, loc="lower left")
axR.set_ylabel("$\\sigma_{resid}$ after affine fit\n(prevalence points)")
axR.grid(alpha=.25)
axA.set_ylabel("Availability\n(% estimable)")
axA.set_ylim(0, 105); axA.grid(alpha=.25)
axA.set_xlabel("Carelessness criterion: minimum number of failed checks ($k$)")
axA.set_xticks([1, 2, 3, 4, 5])
plt.tight_layout(); plt.savefig("strictness_dashboard.png", dpi=160); plt.close()
print("\nwrote strictness_dashboard.png, strictness_dashboard.csv")
