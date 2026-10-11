#!/usr/bin/env python3
"""
Zonation-entropy vs plasticity pipeline  (Gribben GSE202379 / MASLD).

CS framing:
  - Train a 2-class "zone" logistic head on healthy hepatocytes (periportal vs
    pericentral), but the "weights" are just gene-set membership -> robust, no fitting.
  - eta  = P(periportal)  in (0,1)          # the classifier's confidence
  - H    = Shannon entropy of eta            # "how lost is this cell's position"
  - plast= min(hep_score, chol_score)        # expresses BOTH programs -> biphenotypic
  - GROUP BY donor -> GROUP BY stage -> two curves -> 2 stats tests.

--------------------------------------------------------------------------------
STEP 0  (run ONCE in R to export marker counts; Python can't read Seurat .rds):

  library(Seurat)
  obj <- readRDS("GSE202379_SeuratObject_AllCells.rds")   # inner file after un-gzip
  markers <- c("HAL","SDS","IGFBP2","ASS1",
               "GLUL","CYP3A4","SLCO1B3","LGR5","AXIN2",
               "ALB","HNF4A","KRT19","SOX9","EPCAM")
  m <- GetAssayData(obj, slot="counts")[markers, ]        # genes x cells, raw counts
  write.csv(as.matrix(t(m)), "marker_counts.csv")         # cells x genes  (+ rownames=cell id)
--------------------------------------------------------------------------------
"""

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, wilcoxon
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ----------------------------------------------------------------------------- config
META   = "gribben_metadata.csv"     # has Disease.status, cell.annotation, Patient.ID, cell id index
COUNTS = "marker_counts.csv"        # cells x genes, raw counts (from STEP 0)
OUTDIR = "."

# EDIT these to the EXACT strings in Disease.status (case matters), low->high severity:
STAGE_ORDER = ["healthy", "NAFLD", "NASH", "cirrhosis", "end-stage"]

PERIPORTAL  = ["HAL", "SDS", "IGFBP2", "ASS1"]          # eta -> 1  (zone 8)
PERICENTRAL = ["GLUL", "CYP3A4", "SLCO1B3", "LGR5", "AXIN2"]  # eta -> 0  (zone 1)
HEP_PROG    = ["ALB", "HNF4A"]
CHOL_PROG   = ["KRT19", "SOX9", "EPCAM"]

LOW_ENTROPY_Q  = 0.50   # RQ2: keep cells with confident zone call (entropy below this quantile)
HIGH_PLAST_Q   = 0.75   # RQ2: "high-plasticity" = above this quantile
N_BOOT         = 2000   # RQ1 precedence bootstrap over donors
RNG            = np.random.default_rng(0)


# ----------------------------------------------------------------------------- helpers
def sigmoid(x):
    return 1.0 / (1.0 + np.exp(-x))

def binary_entropy(p, eps=1e-9):
    p = np.clip(p, eps, 1 - eps)
    return -(p * np.log2(p) + (1 - p) * np.log2(1 - p))   # in [0,1], max at p=0.5

def zscore(df):
    """z-score each gene (column) across cells. Returns same shape."""
    return (df - df.mean(axis=0)) / (df.std(axis=0) + 1e-9)

def gene_set_score(zexpr, genes):
    """mean z-expression over a gene set -> per-cell score (equal weights)."""
    present = [g for g in genes if g in zexpr.columns]
    return zexpr[present].mean(axis=1)


# ----------------------------------------------------------------------------- load + normalize
def load():
    meta   = pd.read_csv(META, index_col=0)
    counts = pd.read_csv(COUNTS, index_col=0)
    cells  = meta.index.intersection(counts.index)
    meta, counts = meta.loc[cells], counts.loc[cells]

    # Seurat-style log-normalize: counts/total*1e4 -> log1p ; THEN z-score per gene.
    lib = counts.sum(axis=1).replace(0, np.nan)
    logn = np.log1p(counts.div(lib, axis=0) * 1e4).fillna(0.0)
    z = zscore(logn)
    return meta, z


# ----------------------------------------------------------------------------- per-cell features
def per_cell_features(meta, z):
    pp = gene_set_score(z, PERIPORTAL)
    pc = gene_set_score(z, PERICENTRAL)
    eta   = sigmoid(pp - pc)                 # step 4  (logistic head, no fitting)
    ent   = binary_entropy(eta)              # step 5  "lost position"
    hep   = gene_set_score(z, HEP_PROG)
    chol  = gene_set_score(z, CHOL_PROG)
    plast = np.minimum(hep, chol)            # step 6  biphenotypic = high in BOTH

    f = pd.DataFrame({"eta": eta, "entropy": ent, "plasticity": plast}, index=meta.index)
    for c in ["Disease.status", "cell.annotation", "Patient.ID"]:
        f[c] = meta[c].values
    return f


# ----------------------------------------------------------------------------- aggregate
def stage_curve(feat, col):
    """donor-mean -> stage mean +/- SE across donors. Returns df indexed by STAGE_ORDER."""
    donor = feat.groupby(["Disease.status", "Patient.ID"])[col].mean().reset_index()
    g = donor.groupby("Disease.status")[col]
    out = pd.DataFrame({"mean": g.mean(),
                        "se": g.std() / np.sqrt(g.count()),
                        "n_donors": g.count()})
    return out.reindex(STAGE_ORDER)


# ----------------------------------------------------------------------------- RQ1 precedence
def first_cross_stage(donor_means_by_stage):
    """index (>=1) of earliest diseased stage whose mean exceeds healthy baseline; else +inf."""
    base = donor_means_by_stage[STAGE_ORDER[0]]
    for i, s in enumerate(STAGE_ORDER[1:], start=1):
        if donor_means_by_stage.get(s, -np.inf) > base:
            return i
    return np.inf

def rq1_precedence(feat):
    donor = feat.groupby(["Disease.status", "Patient.ID"])[["entropy", "plasticity"]].mean().reset_index()
    donors_by_stage = {s: donor[donor["Disease.status"] == s] for s in STAGE_ORDER}

    earlier = 0
    for _ in range(N_BOOT):
        ent_means, pl_means = {}, {}
        for s in STAGE_ORDER:
            d = donors_by_stage[s]
            if len(d) == 0:
                ent_means[s] = pl_means[s] = -np.inf
                continue
            idx = RNG.integers(0, len(d), len(d))           # resample donors w/ replacement
            ent_means[s] = d["entropy"].values[idx].mean()
            pl_means[s]  = d["plasticity"].values[idx].mean()
        if first_cross_stage(ent_means) < first_cross_stage(pl_means):
            earlier += 1
    p_entropy_first = earlier / N_BOOT

    # point estimate on observed data
    ent_obs = {s: donors_by_stage[s]["entropy"].mean()    for s in STAGE_ORDER}
    pl_obs  = {s: donors_by_stage[s]["plasticity"].mean() for s in STAGE_ORDER}
    return {
        "entropy_first_cross": first_cross_stage(ent_obs),
        "plasticity_first_cross": first_cross_stage(pl_obs),
        "P(entropy crosses earlier)": p_entropy_first,
    }


# ----------------------------------------------------------------------------- RQ2+ zone x stage plasticity
def rq2plus_donor_table(feat):
    """Assign a confident periportal/pericentral label, then per-donor mean plasticity
    per (zone, stage). Returns donor-level long table used by both RQ2+ tests."""
    conf = feat[feat["entropy"] <= feat["entropy"].quantile(LOW_ENTROPY_Q)].copy()
    conf["zone"] = np.where(conf["eta"] >= 0.5, "periportal", "pericentral")
    donor = (conf.groupby(["zone", "Disease.status", "Patient.ID"])["plasticity"]
                 .mean().reset_index())                       # 1 row per (zone, stage, donor)
    return donor


def _zone_first_cross(donor, zone):
    """earliest diseased stage where this zone's donor-mean plasticity > its healthy mean."""
    d = donor[donor["zone"] == zone]
    means = d.groupby("Disease.status")["plasticity"].mean()
    base = means.get(STAGE_ORDER[0], np.inf)
    for i, s in enumerate(STAGE_ORDER[1:], start=1):
        if means.get(s, -np.inf) > base:
            return i
    return np.inf


def rq2plus_zone_plasticity(feat):
    donor = rq2plus_donor_table(feat)
    zones = ["periportal", "pericentral"]

    # --- trajectory: mean +/- SE across donors, per (zone, stage) ---
    traj = {}
    for zn in zones:
        d = donor[donor["zone"] == zn]
        g = d.groupby("Disease.status")["plasticity"]
        traj[zn] = pd.DataFrame({"mean": g.mean(),
                                 "se": g.std() / np.sqrt(g.count()),
                                 "n": g.count()}).reindex(STAGE_ORDER)

    # --- "MORE": at each stage, paired (per-donor) periportal vs pericentral ---
    more = {}
    for s in STAGE_ORDER:
        w = donor[donor["Disease.status"] == s].pivot_table(
            index="Patient.ID", columns="zone", values="plasticity")
        w = w.dropna(subset=zones)                            # donors with both zones present
        if len(w) >= 3:
            try:
                stat, p = wilcoxon(w["periportal"], w["pericentral"])
            except ValueError:
                stat, p = np.nan, np.nan
            higher = "periportal" if w["periportal"].mean() > w["pericentral"].mean() else "pericentral"
            more[s] = {"n_donors": len(w), "p": p, "more_plastic_zone": higher}
        else:
            more[s] = {"n_donors": len(w), "p": np.nan, "more_plastic_zone": "n/a"}

    # --- "EARLIER": which zone crosses its own healthy baseline first (bootstrap) ---
    pp_first = _zone_first_cross(donor, "periportal")
    pc_first = _zone_first_cross(donor, "pericentral")
    pp_ids = donor[donor["zone"] == "periportal"]
    pc_ids = donor[donor["zone"] == "pericentral"]
    earlier_pp = 0
    for _ in range(N_BOOT):
        def boot_first(d):
            means = {}
            for s in STAGE_ORDER:
                v = d[d["Disease.status"] == s]["plasticity"].values
                means[s] = v[RNG.integers(0, len(v), len(v))].mean() if len(v) else -np.inf
            base = means[STAGE_ORDER[0]]
            for i, s in enumerate(STAGE_ORDER[1:], start=1):
                if means[s] > base:
                    return i
            return np.inf
        if boot_first(pp_ids) < boot_first(pc_ids):
            earlier_pp += 1

    earlier = {"periportal_first_cross": pp_first,
               "pericentral_first_cross": pc_first,
               "P(periportal plastic earlier)": earlier_pp / N_BOOT}
    return traj, more, earlier


# ----------------------------------------------------------------------------- plot
def plot_curves(ent_c, pl_c):
    x = range(len(STAGE_ORDER))
    fig, ax1 = plt.subplots(figsize=(7, 4))
    ax1.errorbar(x, ent_c["mean"], yerr=ent_c["se"], marker="o", color="tab:blue", label="zonation entropy")
    ax1.set_ylabel("mean entropy (lost position)", color="tab:blue")
    ax2 = ax1.twinx()
    ax2.errorbar(x, pl_c["mean"], yerr=pl_c["se"], marker="s", color="tab:red", label="plasticity")
    ax2.set_ylabel("mean plasticity (changing type)", color="tab:red")
    ax1.set_xticks(list(x)); ax1.set_xticklabels(STAGE_ORDER, rotation=20)
    ax1.set_title("RQ1: which curve leaves healthy baseline first?")
    fig.tight_layout(); fig.savefig(f"{OUTDIR}/rq1_curves.png", dpi=150)


def plot_zone_trajectory(traj):
    x = range(len(STAGE_ORDER))
    fig, ax = plt.subplots(figsize=(7, 4))
    for zn, color in [("periportal", "tab:green"), ("pericentral", "tab:purple")]:
        ax.errorbar(x, traj[zn]["mean"], yerr=traj[zn]["se"], marker="o", color=color, label=zn)
    ax.set_xticks(list(x)); ax.set_xticklabels(STAGE_ORDER, rotation=20)
    ax.set_ylabel("mean plasticity"); ax.legend()
    ax.set_title("RQ2+: does one zone go plastic earlier/more?")
    fig.tight_layout(); fig.savefig(f"{OUTDIR}/rq2plus_zone_trajectory.png", dpi=150)


# ----------------------------------------------------------------------------- main
def main():
    meta, z = load()
    feat = per_cell_features(meta, z)
    feat.to_csv(f"{OUTDIR}/per_cell_features.csv")        # the one table both RQs use

    hep = feat[feat["cell.annotation"] == "Hepatocytes"]  # RQ1 + RQ2+: hepatocytes (zone only means sth here)

    ent_c = stage_curve(hep, "entropy")
    pl_c  = stage_curve(hep, "plasticity")
    plot_curves(ent_c, pl_c)

    print("=== stage curves (hepatocytes, donor-averaged) ===")
    print("ENTROPY\n", ent_c, "\nPLASTICITY\n", pl_c)
    print("\n=== RQ1 precedence ===");  print(rq1_precedence(hep))

    traj, more, earlier = rq2plus_zone_plasticity(hep)
    plot_zone_trajectory(traj)
    print("\n=== RQ2+ trajectory (periportal) ===\n", traj["periportal"])
    print("\n=== RQ2+ trajectory (pericentral) ===\n", traj["pericentral"])
    print("\n=== RQ2+ MORE: per-stage zone contrast (paired) ===")
    for s in STAGE_ORDER:
        print(f"  {s}: {more[s]}")
    print("\n=== RQ2+ EARLIER ===\n", earlier)


if __name__ == "__main__":
    main()
