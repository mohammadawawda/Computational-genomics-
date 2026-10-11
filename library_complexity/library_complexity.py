import sys, os, glob, re
import numpy as np
from scipy.optimize import brentq
from scipy.stats import chi2, poisson
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

DATA_DIRNAME = "data"


def read_fragments(path):
    """Read one file -> dict {(chr,start,end): occurrences}, summing dups."""
    d = {}
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            p = re.split(r"[,\s]+", line)      # comma, tab, or space delimited
            if len(p) < 4:
                continue
            try:                                # skip malformed lines
                key = (p[0], int(p[1]), int(p[2]))
                cnt = int(p[3])
            except ValueError:
                continue
            d[key] = d.get(key, 0) + cnt
    return d


def combine(paths):
    """Merge several files into length and k arrays """
    merged = {}
    for path in paths:
        for key, cnt in read_fragments(path).items():
            merged[key] = merged.get(key, 0) + cnt
    length = np.array([e - s + 1 for (_, s, e) in merged], dtype=int)
    k = np.array(list(merged.values()), dtype=int)
    return length, k


def estimate_lambda(k):
    """
    Pearson chi-squared test of observed vs predicted counts under the fitted
    positive Poisson. Sparse tail categories are pooled until each expected
    count reaches 5. Returns (stat, dof, p, (js, O, E)) with unmerged arrays.
    """
    M_u = len(k)
    M = k.sum()
    if M_u == 0:
        return np.nan
    r = M_u / M
    if r >= 1.0:
        return np.nan
    g = lambda lam: (1 - np.exp(-lam)) / lam - r
    return brentq(g, 1e-9, 1e6)


def trunc_poisson_pmf(j, lam):
    """P(k=j | k>0) for the positive Poisson, j >= 1."""
    return poisson.pmf(j, lam) / (1 - np.exp(-lam))


def gof_test(k, lam):
    """
    Pearson chi-squared test comparing observed and predicted counts under the
    fitted positive Poisson. Sparse tail categories get pooled until each
     expected count reaches 5. Returns (stat, dof, p, (js, O, E)) with the
    unmerged arrays.
    """
    M_u = len(k)
    js = np.arange(1, k.max() + 1)
    O = np.array([(k == j).sum() for j in js], dtype=float)
    E = M_u * trunc_poisson_pmf(js, lam)

    O_m, E_m = [], []
    co = ce = 0.0
    for o, e in zip(O, E):
        co += o; ce += e
        if ce >= 5:
            O_m.append(co); E_m.append(ce); co = ce = 0.0
    if ce > 0:
        if O_m:
            O_m[-1] += co; E_m[-1] += ce
        else:
            O_m.append(co); E_m.append(ce)
    O_m, E_m = np.array(O_m), np.array(E_m)

    stat = ((O_m - E_m) ** 2 / E_m).sum()
    dof = max(len(O_m) - 1 - 1, 1)
    return stat, dof, chi2.sf(stat, dof), (js, O, E)


# truncated poisson fit and goodness of fit for every sample
def fit_samples(samples):
    for name, (length, k) in samples.items():
        lam = estimate_lambda(k)
        N_hat = len(k) / (1 - np.exp(-lam))
        stat, dof, p, (js, O, E) = gof_test(k, lam)

        print(f"\n=== {name} ===")
        print(f"  M  (total reads)   = {int(k.sum())}")
        print(f"  Mu (unique mols)   = {len(k)}")
        print(f"  lambda_hat         = {lam:.4f}   (same as the MLE)")
        print(f"  N_hat              = {N_hat:.1f}")
        print(f"  chi2 GoF: stat={stat:.2f}, dof={dof}, p={p:.3g}")

        plt.figure()
        plt.bar(js - 0.2, O, width=0.4, label="observed")
        plt.bar(js + 0.2, E, width=0.4, label="predicted")
        plt.xlabel("k (occurrences)"); plt.ylabel("# molecules")
        plt.title(f"{name}  (lambda={lam:.3f}, p={p:.2g})")
        plt.legend(); plt.tight_layout()
        plt.savefig(f"{name}.hist.png", dpi=120); plt.close()


# how the estimate changes with the fragment length
def lambda_by_length(samples, binsize=5, min_mols=30):
    per = {}
    plt.figure()
    for name, (length, k) in samples.items():
        left = (length // binsize) * binsize
        centers, lams, Ns = [], [], []
        for b in np.unique(left):
            kk = k[left == b]
            if len(kk) < min_mols:
                continue
            lam = estimate_lambda(kk)
            if np.isnan(lam):
                continue
            centers.append(b + binsize / 2)
            lams.append(lam)
            Ns.append(len(kk) / (1 - np.exp(-lam)))
        per[name] = dict(length=np.array(centers),
                         lam=np.array(lams), N=np.array(Ns))
        plt.plot(centers, lams, marker="o", label=name)
    plt.xlabel("fragment length (bp)"); plt.ylabel("estimated lambda")
    plt.title("lambda vs fragment length"); plt.legend(); plt.tight_layout()
    plt.savefig("lambda_vs_length.png", dpi=120); plt.close()

    plt.figure()
    for name, d in per.items():
        plt.plot(d["length"], d["N"], marker="o", label=name)
    plt.xlabel("fragment length (bp)")
    plt.ylabel("estimated original molecules  N(L)")
    plt.title("Reconstructed original fragment-length distribution")
    plt.legend(); plt.tight_layout()
    plt.savefig("original_length_dist.png", dpi=120); plt.close()


def discover_samples():
    """Find the data folder next to this script and group files by sample."""
    here = os.path.dirname(os.path.abspath(__file__))
    data_dir = os.path.join(here, DATA_DIRNAME)
    if not os.path.isdir(data_dir):
        return None
    files = sorted(f for f in glob.glob(os.path.join(data_dir, "*"))
                   if os.path.isfile(f))
    groups = {}
    for f in files:
        base = os.path.basename(f)

        name = os.path.splitext(base)[0]
        groups.setdefault(name, []).append(f)
    return groups


if __name__ == "__main__":

    groups = discover_samples()
    if not groups:
        sys.exit(f'No files. Put data in a "{DATA_DIRNAME}" folder next to '
                 f'this script, or pass file paths as arguments.')

    print("Samples found:")
    for name, paths in groups.items():
        print(f"  {name}: {[os.path.basename(p) for p in paths]}")

    samples = {name: combine(paths) for name, paths in groups.items()}
    fit_samples(samples)
    lambda_by_length(samples)
    print("\nPlots: <sample>.hist.png, lambda_vs_length.png, original_length_dist.png")
