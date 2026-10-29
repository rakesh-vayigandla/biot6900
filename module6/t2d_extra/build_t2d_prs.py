"""Module 6 extra (not graded): the same portability workflow on a REAL type 2 diabetes score.

Why: my Part I disease is T2D, but the course genotype file contains only the 22 Alzheimer's
score SNPs, so a T2D score cannot be run on it. Here I extract the genotypes a T2D score
actually needs, straight from the 1000 Genomes VCFs, for the same 2,336 people.

Inputs  : PGS Catalog PGS003733 (EUR-derived) and PGS003735 (AFR-derived), DIAMANTE T2D;
          1000 Genomes phase 3 GRCh37 VCFs (remote); ../data/meta.csv for the sample list.
Outputs : t2d_prs_by_ancestry.csv, figures/t2d_prs_by_ancestry.png

There is no T2D phenotype for 1000 Genomes, so this reports the PRS DISTRIBUTION only.
No AUC or calibration is computed, and none is implied.

Reproduce: bcftools (conda install -c bioconda bcftools), then
  python build_t2d_prs.py          # re-extraction step is in raw/ (see README note)
"""
import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt

SCORES = {"EUR-derived (PGS003733)": "raw/PGS003733_hmPOS_GRCh37.txt.gz",
          "AFR-derived (PGS003735)": "raw/PGS003735_hmPOS_GRCh37.txt.gz"}
POPS = ["EUR", "AFR", "EAS", "SAS", "AMR"]
COL = {"EUR": "#1E2028", "AFR": "#C8102E", "EAS": "#6E7178", "SAS": "#A4804A", "AMR": "#9AA0A6"}

samples = [l.strip() for l in open("raw/vcf_samples.txt")]
rows = []
for line in open("raw/geno_raw.tsv"):                       # from bcftools, see module README
    f = line.rstrip("\n").split("\t")
    gt = [g.count("1") if g[0] in "01" and g[-1] in "01" else np.nan for g in f[5:]]
    rows.append((f[0], int(f[1]), f[3], f[4], gt))
G = pd.DataFrame({(r[0], r[1]): r[4] for r in rows}, index=samples)
alleles = {(r[0], r[1]): (r[2], r[3]) for r in rows}        # (REF, ALT) per position
meta = pd.read_csv("../data/meta.csv")
assert list(meta.sample_id) == samples, "sample order must match meta.csv"

prs = {}
for name, path in SCORES.items():
    s = pd.read_csv(path, sep="\t", comment="#").dropna(subset=["hm_chr", "hm_pos"])
    s["hm_chr"] = s.hm_chr.astype(str); s["hm_pos"] = s.hm_pos.astype(int)
    used = flip = drop = absent = 0
    total = np.zeros(len(G))
    for _, r in s.iterrows():                               # keep / flip / drop, as in the graded notebook
        k = (r.hm_chr, r.hm_pos)
        if k not in G.columns:
            absent += 1; continue
        ref, alt = alleles[k]
        ea = str(r.effect_allele).upper()
        d = G[k].to_numpy(float)
        if ea == alt:
            used += 1
        elif ea == ref:
            d = 2 - d; flip += 1; used += 1
        else:
            drop += 1; continue
        total += np.nan_to_num(d, nan=np.nanmean(d)) * r.effect_weight
    prs[name] = (total - total.mean()) / total.std()
    print(f"{name}: {len(s)} variants | used {used} (flipped {flip}) | allele-mismatch drops {drop} | not in panel {absent}")

out = meta[["sample_id", "ancestry"]].copy()
for k, v in prs.items():
    out[k.split()[0] + "_prs"] = np.round(v, 4)
out.to_csv("t2d_prs_by_ancestry.csv", index=False)

eur = prs["EUR-derived (PGS003733)"]
cut = np.quantile(eur, 0.90)
print(f"\nwhole-cohort top-10% cutoff on the EUR-derived score: {cut:.2f}")
for p in POPS:
    m = (meta.ancestry == p).to_numpy()
    print(f"  {p}: mean {eur[m].mean():+.2f} (SD {eur[m].std():.2f}) | above cutoff {100*(eur[m] > cut).mean():.1f}%")
print(f"\ncorrelation between the EUR- and AFR-derived scores: {np.corrcoef(*prs.values())[0,1]:.3f}")

import os; os.makedirs("figures", exist_ok=True)
fig, axes = plt.subplots(1, 2, figsize=(10, 3.8), sharex=True, sharey=True)
for ax, (name, z) in zip(axes, prs.items()):
    for p in POPS:
        ax.hist(z[(meta.ancestry == p).to_numpy()], bins=28, density=True, alpha=0.4, color=COL[p], label=p)
    ax.set_title(name, fontsize=10); ax.set_xlabel("PRS (standardized)")
    for sp in ("top", "right"): ax.spines[sp].set_visible(False)
axes[0].set_ylabel("Density"); axes[0].axvline(cut, color="k", ls="--", lw=1)
axes[0].annotate("whole-cohort\ntop 10%", (cut, 0.42), fontsize=7, ha="left")
axes[1].legend(frameon=False, fontsize=8)
fig.suptitle("Real type 2 diabetes scores on 1000 Genomes (no phenotype: distribution only)", fontsize=10.5)
plt.tight_layout(); plt.savefig("figures/t2d_prs_by_ancestry.png", dpi=200)
print("\nwrote t2d_prs_by_ancestry.csv and figures/t2d_prs_by_ancestry.png")
