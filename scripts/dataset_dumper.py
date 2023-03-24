"""
Dump dataset to LOGAN format and perform some statistics.
"""
#%%
from dataclasses import dataclass
import json
import gzip

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

# %%

# ecoli = Dataset.from_json(*ECOLI_ARGS)
# dump_logan(ecoli, OUTDIR + "ecoli.tsv.gz")

# %%

datasets = [
    # ("simulated85", SIM85_ARGS),
    # ("simulated0", SIM0_ARGS),
    # ("simulated1", SIM1_ARGS),
    ("ecoli", ECOLI_ARGS),
    ("ecoli100", ECOLI100_ARGS),
    ("elegans", ELEGANS_ARGS),
]
for name, args in datasets:
    ds = Dataset.from_json(*args)
    dump_logan(ds, OUTDIR + f"{name}.tsv.gz")

# %%
