#%matplotlib notebook
"""
Dump dataset to LOGAN format and perform some statistics.
"""
from pathlib import Path
from argparse import ArgumentParser
from dataclasses import dataclass
import json
import gzip
from pprint import pprint

import numpy as np
import pandas as pd

def load_json(p):
    p = str(p)
    if p.endswith(".gz"):
        with gzip.open(p, "r") as gf:
            data = json.load(gf)
    else:
        with open(p) as f:
            data = json.load(f)
    return data

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
        "seqcount": len(ds.sequences),
        "seqlen_avg": np.average(seq_l + seq_r),
        "seqlen_avg_l": np.average(seq_l),
        "seqlen_avg_r": np.average(seq_r),
        "seqlen_std": np.std(seq_l + seq_r),
        "seqlen_std_l": np.std(seq_l),
        "seqlen_std_r": np.std(seq_r),
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


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("--seedlength", type=int, default="17")
    parser.add_argument("cmps", type=Path)
    parser.add_argument("seqs", type=Path)
    args = parser.parse_args()
    ds = Dataset.from_json(args.seqs, args.cmps)
    stats, *_ = get_dataset_stats("dataset", ds, seedlen=args.seedlength)
    pprint(stats)
