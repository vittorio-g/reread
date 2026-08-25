# plot_injection_curve.py — estimated vs true careless prevalence, by injection.
# Solid = true prevalence of the constructed sample; dotted = prevalence estimated by CaReReRe.
# Faceted by dataset family so no panel carries more than three series (an 8-series categorical
# palette cannot pass CVD separation; faceting is the recommended fix). Each series also carries a
# distinct marker and a direct end-label, the secondary encoding the palette's dE band requires.
import csv, collections
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
plt.rcParams.update({"font.family": "serif", "font.size": 10.5})

rows = list(csv.DictReader(open("injection_curve.csv")))
data = collections.OrderedDict()
for r in rows:
    data.setdefault(r["dataset"], []).append(
        (100*float(r["true_prev"]), 100*float(r["pi_mean"]), 100*float(r["pi_sd"]), int(r["n"])))

PANELS = [("TISP (68-country study, by country)",
           ["TISP Germany", "TISP Poland", "TISP Australia"]),
          ("SMARVUS (multi-site study, by country)",
           ["smarvus England", "smarvus Egypt", "smarvus Canada"]),
          ("Other datasets",
           ["warning control", "Kay S1", "TUMI"])]
COL = ["#B35806", "#2166AC", "#762A83"]      # validated: CVD dE 7.5 -> markers + labels required
MK  = ["o", "s", "^"]

fig, axes = plt.subplots(1, 3, figsize=(13.6, 4.8), sharex=True, sharey=True)
for ax, (title, names) in zip(axes, PANELS):
    # the constructed true prevalence IS the x coordinate, identical for every series here,
    # so it is drawn once as a single reference line rather than as three overlapping ones
    ax.plot([0, 42], [0, 42], color="#555555", lw=2.0, ls="-", zorder=2)
    ax.text(39.5, 41.5, "true", color="#555555", fontsize=8.5, ha="right", rotation=40)
    for i, nm in enumerate(names):
        d = data.get(nm)
        if not d: continue
        x = [p[0] for p in d]; est = [p[1] for p in d]
        ax.plot(x, est, color=COL[i], lw=2.0, ls=":", marker=MK[i], ms=5,
                zorder=3, label=nm.split(" ", 1)[-1])                              # estimated
        if est:
            ax.annotate(nm.split(" ", 1)[-1], xy=(x[-1], est[-1]), xytext=(3, 0),
                        textcoords="offset points", color=COL[i], fontsize=8.5,
                        va="center", fontweight="bold")
    ax.set_title(title, fontsize=10.5)
    ax.set_xlabel("True careless rate (%)")
    ax.grid(alpha=.25); ax.set_xlim(0, 46); ax.set_ylim(0, 46)
axes[0].set_ylabel("Prevalence (%)")
from matplotlib.lines import Line2D
axes[0].legend(handles=[Line2D([0],[0], color="#555555", lw=2, ls="-", label="true prevalence (constructed)"),
                        Line2D([0],[0], color="#555555", lw=2, ls=":", marker="o", ms=5, label=r"estimated $\hat\pi$")],
               loc="upper left", fontsize=8.5, framealpha=.92)
fig.suptitle("Estimated versus true careless prevalence: real respondents injected at controlled rates "
             "(10 resamples per point)", fontsize=11.5, y=1.00)
plt.tight_layout()
plt.savefig("injection_curve.png", dpi=160, bbox_inches="tight"); plt.close()
print("wrote injection_curve.png")
