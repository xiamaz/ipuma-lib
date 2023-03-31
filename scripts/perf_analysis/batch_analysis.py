from pathlib import Path

from argparse import ArgumentParser

import pandas as pd
import ipulog_parser

from importlib import reload

reload(ipulog_parser)



parser = ArgumentParser()
parser.add_argument("logdir", type=Path)
parser.add_argument("output_path", type=Path)

args = parser.parse_args()

logdir = args.logdir

outdir = Path("./output/elba_scale/plots")
outdir.mkdir(exist_ok=True)

logdir = Path("./output/elba_scale/elba_scale")
logfiles = [
    ipulog_parser.load_logfile(log) for log in logdir.iterdir()
]

all_df = pd.DataFrame.from_records(logfiles)
all_df["Number IPU"] = all_df["name"].apply(lambda n: int(n.split("ipu")[-1]))
all_df = all_df.sort_values(by="Number IPU", axis=0).reset_index(drop=True)

all_df[["Number IPU", "total_execution_time_ms"]]

import proplot as pplt


def plot_per_ipu_scaling(all_df):
    plot_df=all_df.loc[all_df["Number IPU"] != 64].copy()
    single_perf = all_df.iloc[0]["total_execution_time_ms"]
    theoretical_scaling = [
        single_perf / i for i in range(1, 17)
    ]
    fig = pplt.figure()
    ax = fig.add_subplot()
    ax.format(yscale="log", xlocator=pplt.Locator("fixed", plot_df["Number IPU"]))
    ax.line(
        x="Number IPU",
        y="total_execution_time_ms",
        data=plot_df,
        label="IPU Scaling",
    )
    ax.line(
        y=theoretical_scaling,
        x=range(1, 17),
        label="Theoretical linear scaling",
    )

    fig.suptitle("Scaling 1-16 IPUs on the E coli dataset ELBA")
    ax.legend()
    outpath = outdir / "IPU_scaling_ecoli.png"
    fig.savefig(str(outpath), dpi=300)


def plot_multicomparison_histogram(all_df):
    histogram = pd.DataFrame(all_df.iloc[0]["histogram"].values(), index=all_df.iloc[0]["histogram"].keys(), columns=["Number of comparisons"]).astype(int)
    histogram = histogram.loc[histogram["Number of comparisons"] > 0]
    fig = pplt.figure()
    fig.suptitle("Comparison merging by shared input sequences")
    ax = fig.add_subplot()
    ax.bar(histogram)
    ax.format(xlocator="auto", ylabel="Number of items", xlabel="Merged comparison count")
    fig.savefig(outdir / "IPU_histogram_multicomparison.png", dpi=300)


def plot_gcups_per_batch(all_df):
    all_dfs = []
    for ipucount, logs in zip(all_df["Number IPU"], all_df["entries"]):
        df = pd.DataFrame.from_records([l["data"] for l in logs if l["type"] == "JOBLOG"])
        df["Number IPU"] = ipucount
        df["Number IPU Label"] = f"{ipucount}"
        all_dfs.append(df)


    joblog_df = pd.concat(all_dfs)
    plot_df = joblog_df.pivot(columns="Number IPU Label", values="gcups_outer")
    plot_df = plot_df[map(str, sorted([int(c) for c in plot_df.columns.tolist()]))]
    fig = pplt.figure()
    ax = fig.add_subplot()

    ax.boxplot(plot_df)
    ax.format(xlabel="Number of IPUs", ylabel="Host measured GCUPS per batch")
    fig.savefig(outdir / "IPU_gcups_outer_num_ipu.png", dpi=300)


plot_per_ipu_scaling(all_df)
plot_multicomparison_histogram(all_df)
plot_gcups_per_batch(all_df)
