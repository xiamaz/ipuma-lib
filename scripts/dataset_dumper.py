#%matplotlib notebook
"""
Dump dataset to LOGAN format and perform some statistics.
"""
#%%
from dataclasses import dataclass
import json
import gzip

import numpy as np
import pandas as pd
import proplot as pplt

def load_json(p):
    with open(p) as f:
        return json.load(f)

def load_type(p, t):
    data = []
    with open(p) as f:
        for l in f:
            data.append(t(l))
    return data

@dataclass(frozen=True)
class Comparison:
    originalComparisonIndex: int
    indexA: int
    sizeA: int
    indexB: int
    sizeB: int
    complexity: int
    seeds: list = None

@dataclass(frozen=True)
class MultiComparison:
    comparisonCount: int = 0
    comparisons: list = None
    totalSeqSize: int = 0

    @classmethod
    def from_json(cls, data):
        return cls(
            comparisonCount=data["comparisonCount"],
            comparisons=list(map(lambda d: Comparison(**d), data["comparisons"])),
            totalSeqSize=data["totalSeqSize"],
        )

class Dataset:
    def __init__(self, sequences, cmps):
        self.sequences = sequences
        self.cmps = cmps

    @classmethod
    def from_json(cls, seqfile, cmpfile):
        return cls(load_json(seqfile), list(map(MultiComparison.from_json, load_json(cmpfile))))

def dump_logan(dataset: Dataset, output: str):
    print("Writint LOGAN format output to", output)
    with gzip.open(output, "wb", compresslevel=6) as f:
        for mcmp in dataset.cmps:
            for cmp in mcmp.comparisons:
                for seedA, seedB in cmp.seeds:
                    if seedA >= 0 and seedB >= 0:
                        outline = f"{dataset.sequences[cmp.indexA]}\t{seedA}\t{dataset.sequences[cmp.indexB]}\t{seedB}\tn\n"
                        f.write(outline.encode())

ECOLI_ARGS = (
    "/global/D1/projects/ipumer/inputs_ab/seqs_ecoli_multi_single.json",
    "/global/D1/projects/ipumer/inputs_ab/cmps_ecoli_multi_single.json",
)

ECOLI100_ARGS = (
    "/global/D1/projects/ipumer/inputs_ab/seqs_elba100_multi.json",
    "/global/D1/projects/ipumer/inputs_ab/cmps_elba100_multi.json",
)

ELEGANS_ARGS = (
    "/global/D1/projects/ipumer/inputs_ab/seqs_celegans_multi.json",
    "/global/D1/projects/ipumer/inputs_ab/cmps_celegans_multi.json",
)

SIM85_ARGS = (
    "/global/D1/projects/ipumer/datasets_giulia/seqs_simulated85_multi.json",
    "/global/D1/projects/ipumer/datasets_giulia/cmps_simulated85_multi.json",
)

SIM0_ARGS = (
    "/global/D1/projects/ipumer/datasets_giulia/seqs_simulated0_multi.json",
    "/global/D1/projects/ipumer/datasets_giulia/cmps_simulated0_multi.json",
)

SIM1_ARGS = (
    "/global/D1/projects/ipumer/datasets_giulia/seqs_simulated1_multi.json",
    "/global/D1/projects/ipumer/datasets_giulia/cmps_simulated1_multi.json",
)

OUTDIR = "/global/D1/projects/ipumer/datasets_giulia/"

def get_dataset_stats(name, ds, seedlen = 17):
    cmpcount = 0
    seq_r = []
    seq_l = []
    comparison_complexity = []
    for multicomparison in ds.cmps:
        for comparison in multicomparison.comparisons:
            for pH, pV in comparison.seeds:
                if pH >= 0 and pV >= 0:
                    cmpcount += 1
                    H = len(ds.sequences[comparison.indexA])
                    V = len(ds.sequences[comparison.indexB])
                    lH = pH
                    lV = pV
                    rH = H - (pH + seedlen)
                    rV = V - (pV + seedlen)
                    if lH > 0 and lV > 0:
                        seq_l.append(lH)
                        seq_l.append(lV)
                        comparison_complexity.append(lH * lV)
                    if rH > 0 and rV > 0:
                        seq_r.append(rH)
                        seq_r.append(rV)
                        comparison_complexity.append(rH * rV)

    stats = {
        "name": name,
        "cmpcount": cmpcount,
        "seqlen_avg": np.average(seq_l + seq_r),
        "seqlen_avg_l": np.average(seq_l),
        "seqlen_avg_r": np.average(seq_r),
        "seqlen_max_l": np.max(seq_l),
        "seqlen_max_r": np.max(seq_r),
        "seqlen_min_l": np.min(seq_l),
        "seqlen_min_r": np.min(seq_r),
        "seqlen_p90_l": np.percentile(seq_l, 90),
        "seqlen_p90_r": np.percentile(seq_r, 90),
        "seqlen_p10_l": np.percentile(seq_l, 10),
        "seqlen_p10_r": np.percentile(seq_r, 10),
        "complexity": np.mean(comparison_complexity),
    }
    return stats, seq_r, seq_l, comparison_complexity


datasets = [
    ("simulated85", SIM85_ARGS),
    # ("simulated0", SIM0_ARGS),
    # ("simulated1", SIM1_ARGS),
    ("ecoli", ECOLI_ARGS),
    ("ecoli100", ECOLI100_ARGS),
    ("elegans", ELEGANS_ARGS),
]
#%%
dsstats = []
dsvalues = {}
for name, args in datasets:
    ds = Dataset.from_json(*args)
    dsvalues[name] = {}
    stats, dsvalues[name]["sr"], dsvalues[name]["sl"], dsvalues[name]["c"] = get_dataset_stats(name, ds)
    dsstats.append(stats)

#%%
stats = pd.DataFrame.from_records(dsstats)
stats

dfstats_path = "/global/D1/projects/ipumer/dataset_stats.pkl"
stats.to_pickle(dfstats_path)

stats.columns
stats_table = stats[["name", "cmpcount", "seqlen_avg", "seqlen_p10_l", "seqlen_avg_l", "seqlen_p90_l", "seqlen_p10_r", "seqlen_avg_r", "seqlen_p90_r", "complexity"]]
stats_table = stats_table.set_axis(["Name", "Cmp Count", "Seqlen Avg", "Seqlen P10 L", "Seqlen Avg L", "Seqlen P90 L", "Seqlen P10 R", "Seqlen Avg R", "Seqlen P90 R", "Complexity"], axis=1, inplace=False)
def fmt(v):
    if isinstance(v, str):
        return v.replace("_", "\\_")
    else:
        return f"${v:,.0f}$".replace(",", "\\,")
l = stats_table.style.format(fmt, precision=0).to_latex(hrules=True)
print(l)
# CMP COUNT, avg sequence length + max + min R L, avg a * b seed complexity

# plot r and l distribution of comparisons
#%%
import proplot as pplt
fig = pplt.figure(share=False, refaspect=2)
axes = fig.subplots([[1, 2], [3, 4]])

# ax = axes[0]
# ax.hist(testval['sr'][::2], bins=100, label="right", a=.2)
# ax.hist(testval['sr'][1::2], bins=100, label="left", a=.2)
#
# ax = axes[1]
# ax.hist(testval['sl'][::2], bins=100, label="right", a=.2)
# ax.hist(testval['sl'][1::2], bins=100, label="left", a=.2)

cycle = pplt.Colormap("tab10", discrete=True)

patches = []
lower = None
upper = None
for i, name in enumerate(("ecoli", "ecoli100", "elegans")):
    testval = dsvalues[name]
    dslower = min(testval['c'])
    dsupper = max(testval['c'])
    if lower is None:
        lower = dslower
    if upper is None:
        upper = dsupper
    lower = min(dslower, lower)
    upper = max(dsupper, upper)

bins = np.logspace(np.log10(lower), np.log10(upper), 100)

for i, name in enumerate(("ecoli", "ecoli100", "elegans")):
    testval = dsvalues[name]
    if name == "elegans":
        name = "celegans"
    ax = axes[0]
    ax.hist([l + r + 17 for l, r in zip(testval['sl'], testval['sr'])], bins=100, label=name, histtype="step", color=cycle(i), density=True)
    ax.format(title="Total Sequence Length", yticklabels=[])

    ax = axes[1]
    ax.hist(testval['sr'] + testval['sl'], bins=100, label=name, color=cycle(i), histtype="step", density=True)
    ax.format(title="Aligned Sequence Length", yticklabels=[])

    ax = axes[2]
    ratios = [(l + 8.5) / (l + r + 17) for l, r in zip(testval['sl'], testval['sr'])]
    ax.hist(ratios, bins=100, label=name, color=cycle(i), histtype="step", density=True)
    ax.format(title="Relative Seed Position", yticklabels=[])
    lax = ax

    ax = axes[3]
    density, hbins = np.histogram(testval['c'], bins)
    a = ax.stairs(density / density.sum(), hbins, label=name, color=cycle(i))
    ax.format(xscale="log", title="Alignment complexity", yticklabels=[])
    patches.append(a)

lax.legend(patches, loc='b', lw=2, frameon=False)
fig.savefig("/home/zhaom/ipuma-lib/output/dataset_histogram.pdf")


# %%
