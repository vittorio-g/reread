# Who is responsible for the pi-hat collapse? Three panels, one message each:
#  (A) pi-hat vs truth, per index      -> which estimator collapses
#  (B) AUC vs truth, per index         -> is the SIGNAL lost, or only the ESTIMATOR?
#  (C) absolute anchor (raw rc <= 1.5) -> a composition-free quantity that stays monotone
import pandas as pd, numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
plt.rcParams.update({"font.family": "serif", "font.size": 10.5})

d = pd.read_csv("calib_decomp.csv").sort_values("true")
x = d["true"].values * 100

# paper's established categorical order; navy<->rust sits in the 6-8 CVD band, so every
# series also carries a distinct marker + a direct end-label (secondary encoding).
SER = [("eta",          "#8f3535", "o", "ensemble $\\eta$"),
       ("rc",           "#33475f", "s", "$rc$"),
       ("longstring",   "#b08a3a", "^", "LongString"),
       ("person_total", "#2f6b3a", "D", "Person-Total")]

fig, ax = plt.subplots(1, 3, figsize=(13.2, 4.5))

# ---- (A) pi-hat per index -------------------------------------------------
ax[0].plot([0, 47], [0, 47], ls="--", lw=1.2, color="#9a9a9a", zorder=1, label="identity")
for k, c, mk, lab in SER:
    ax[0].plot(x, d["pi_" + k] * 100, marker=mk, ms=5, lw=2, color=c, label=lab, zorder=3)
ax[0].set_xlabel("True careless rate (%)"); ax[0].set_ylabel(r"Estimated $\hat\pi$ (%)")
ax[0].set_title("(a) Only the ENSEMBLE's estimate collapses;\n"
                r"$rc$ and Person-Total stay monotone", fontsize=10.5)
ax[0].grid(alpha=.25); ax[0].legend(fontsize=8.2, loc="upper left", framealpha=.92)
ax[0].set_xlim(0, 47)

# ---- (B) AUC per index ----------------------------------------------------
for k, c, mk, lab in SER:
    ax[1].plot(x, d["auc_" + k], marker=mk, ms=5, lw=2, color=c, label=lab, zorder=3)
ax[1].axhline(.5, ls=":", lw=1, color="#9a9a9a")
ax[1].set_xlabel("True careless rate (%)"); ax[1].set_ylabel("AUC vs true label")
ax[1].set_title("(b) …while every index keeps its ranking signal\n(the failure is the estimator, not the indices)", fontsize=10.5)
ax[1].grid(alpha=.25); ax[1].legend(fontsize=8.2, loc="lower left", framealpha=.92)
ax[1].set_xlim(0, 47); ax[1].set_ylim(.35, 1.02)

# ---- (C) a computable regime test -----------------------------------------
diff = (d["pi_rc"] - d["pi_eta"]).values * 100
ax[2].axhspan(-12, 0, color="#2f6b3a", alpha=.07)
ax[2].axhspan(0, 28, color="#8f3535", alpha=.07)
ax[2].axhline(0, color="#555", lw=1.2, ls="--")
ax[2].plot(x, diff, marker="o", ms=5.5, lw=2.2, color="#33475f", zorder=3)
ax[2].annotate(r"$\hat\pi_{rc} < \hat\pi_{\eta}$: identified regime", xy=(2, -9.2),
               fontsize=8.6, color="#2f6b3a")
ax[2].annotate(r"$\hat\pi_{rc} > \hat\pi_{\eta}$: collapsed regime" "\n(true rate is high)",
               xy=(2, 19), fontsize=8.6, color="#8f3535")
ax[2].set_xlabel("True careless rate (%)")
ax[2].set_ylabel(r"$\hat\pi_{rc} - \hat\pi_{\eta}$ (percentage points)")
ax[2].set_title("(c) The disagreement between the two flips sign\n"
                "at ~25%: a label-free regime test", fontsize=10.5)
ax[2].grid(alpha=.25); ax[2].set_xlim(0, 47); ax[2].set_ylim(-12, 28)

plt.tight_layout()
OUT = "calib_decomp.png"
plt.savefig(OUT, dpi=160); plt.close()
print("wrote " + OUT)
