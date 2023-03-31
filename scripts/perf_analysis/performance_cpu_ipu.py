from pathlib import Path
import pandas as pd
import ipulog_parser
import re

import proplot as pplt

datasets = {
    "cpu": Path("./output/performance_comparison/cpusw_benchmark"),
    "ipu": Path("./output/performance_comparison/ipusw_benchmark_final"),
}

CPUNAMES = {
    "g002": "AMD EPYC 7763",
    "h001": "Intel Xeon 8360Y",
    "ipu-pod64-server1": "Intel Xeon 8168",
}

def format_devicename(devicetype, hostname):
    if devicetype == "cpu":
        return CPUNAMES[hostname]
    elif devicetype == "ipu":
        return "Graphcore IPU"

dataset_stats = {
    logfile.stem: ipulog_parser.load_logfile(logfile)["entries"][-1]["data"]
    for logfile in Path("./output/dataset_stats").glob("*.log")
}

dataset_stats["celegans"] = dataset_stats["elegans"]

def load_gpu_stats(logfile):
    results = []
    with open(logfile) as f:
        dsname = ""
        xvalue = ""
        exectime_inner = 0
        exectime_outer = 0
        exectime_total = 0
        newds = True
        for line in f:
            if line.startswith("#"):
                if dsname:
                    results.append({
                        "device": "gpu",
                        "devicename": "NVIDIA A100",
                        "host": "g002",
                        "dataset": dsname,
                        "xdrop": xvalue,
                        "time_ms": exectime_inner,
                        "time_ms_outer": exectime_outer,
                        "time_ms_total": exectime_total,
                    })
                dsname = line.strip()[1:]
                newds = True
            elif line.startswith("X"):
                if not newds:
                    results.append({
                        "device": "gpu",
                        "devicename": "NVIDIA A100",
                        "host": "g002",
                        "dataset": dsname,
                        "xdrop": xvalue,
                        "time_ms": exectime_inner,
                        "time_ms_outer": exectime_outer,
                        "time_ms_total": exectime_total,
                    })
                else:
                    newds = False
                xvalue = int(line.strip()[2:])
            elif "with data" in line:
                exectime_outer = float(re.search(r"\d+.\d+", line).group(0)) * 1e3
            elif "no data" in line:
                exectime_inner = float(re.search(r"\d+.\d+", line).group(0)) * 1e3
            elif "host execution" in line:
                exectime_total = float(re.search(r"\d+.\d+", line).group(0)) * 1e3
        results.append({
            "device": "gpu",
            "devicename": "NVIDIA A100",
            "host": "g002",
            "dataset": dsname,
            "xdrop": xvalue,
            "time_ms": exectime_inner,
            "time_ms_outer": exectime_outer,
            "time_ms_total": exectime_total,
        })

    return results

gpu_data_giulia = pd.read_csv("./output/performance_comparison/gpu_benchmark/ELBA - Logan - GPU (Perlmutter) - LOGAN-Only.csv").dropna(how="all").drop(columns=["K", "Datasetname"])
gpu_data_giulia["Device"] = gpu_data_giulia["Device"] * 1000
gpu_data_giulia["Host + Device"] = gpu_data_giulia["Host + Device"] * 1000
gpu_data_giulia.columns = ["dataset", "time_ms_total", "time_ms", "devices", "xdrop"]
gpu_data_giulia["devicename"] = "NVIDIA A100"
gpu_data_giulia["algo"] = "logan"
gpu_data_giulia["device"] = "gpu"
gpu_data_giulia["host"] = "perlmutter"
gpu_data_giulia["name"] = gpu_data_giulia.apply(lambda s: f"{s['dataset']}_logan_{s['xdrop']}_gpu{int(s['devices'])}", axis=1)

# gpu_data = pd.DataFrame.from_records(load_gpu_stats("./output/performance_comparison/gpu_benchmark/luks_logan_data.txt"))
all_data = pd.DataFrame.from_records([
    {
        "device": devicetype,
        "devicename": format_devicename(devicetype, hostname.name),
        "host": hostname.name,
        "name": logfile.stem,
        "data": ipulog_parser.load_logfile(logfile),
        "dataset": logfile.stem.split("_")[0],
        "algo": logfile.stem.split("_")[1],
        "xdrop": int(re.search(r"_x(\d+)_", logfile.stem).group(1)), # type: ignore
    }
    for devicetype, dspath in datasets.items()
    for hostname in dspath.iterdir() if hostname.is_dir()
    for logfile in hostname.glob("*.log")
])

def gather_inner_time(s: pd.Series):
    if s["device"] == "cpu" or s["device"] == "gpu":
        return s["time_ms"]

    cycles = sum(e["data"]["cycles_inner"] for e in s["data"]["entries"] if e["type"] == "JOBLOG")
    return cycles / (1.85 * 1e6)

def gather_inner_gcups(s: pd.Series):
    if s["device"] == "cpu":
        if s["dataset"] in ("simulated0", "simulated1"):
            return s["gcups"]
        gcells_ds = dataset_stats[s["dataset"]]["gCellsNoSeed"]
        gcups_computed = gcells_ds / (s["time_ms"] / 1e3)
        return gcups_computed
    if s["device"] == "gpu":
        gcells_ds = dataset_stats[s["dataset"]]["gCellsNoSeed"]
        gcups_computed = gcells_ds / (s["time_ms"] / 1e3)
        return gcups_computed

    cycles = sum(e["data"]["cycles_inner"] for e in s["data"]["entries"] if e["type"] == "JOBLOG")
    if s["dataset"] not in dataset_stats:
        gcells_ds = sum(e["data"]["cell_count"] for e in s["data"]["entries"] if e["type"] == "JOBLOG") / 1e9
    else:
        gcells_ds = dataset_stats[s["dataset"]]["gCellsNoSeed"]
    # cells = sum(e["data"]["cell_count"] for e in s["data"]["entries"] if e["type"] == "JOBLOG")
    gcups_inner_total = (gcells_ds)/(cycles / (1.85*1e9))
    return gcups_inner_total

gpu_data_giulia["gcups_inner"] = gpu_data_giulia.apply(gather_inner_gcups, axis=1)
gpu_data_giulia["time_ms_inner"] = gpu_data_giulia.apply(gather_inner_time, axis=1)
# gpu_data["gcups"] = gpu_data.apply(lambda s: dataset_stats[s["dataset"]]["gCellsNoSeed"] / (s["time_ms_total"] / 1e3), axis=1)

all_data["dataset"].replace({"elegans": "celegans"}, inplace=True)
all_data["time_ms"] = all_data["data"].apply(lambda s: s["entries"][-1]["data"].get("time_ms", None))
# all_data["gcups"] = all_data["data"].apply(lambda s: s["entries"][-1]["data"].get("gcups", None))
all_data["time_ms_inner"] = all_data.apply(gather_inner_time, axis=1)
all_data["gcups_inner"] = all_data.apply(gather_inner_gcups, axis=1)
all_data["threads"] = all_data["data"].apply(lambda s: s["entries"][0]["data"].get("algoconfig", {"threads": 1472})["threads"])
all_data["devices"] = all_data["data"].apply(lambda s: s["entries"][0]["data"].get("numDevices", 1))

not_valid = all_data.loc[all_data["time_ms"].isna()]
if not_valid.shape[0] > 0:
    print("Missing results for following experiments:")
    print(not_valid["name"])


all_valid = all_data.loc[all_data["time_ms"].notna()].copy()
all_valid = pd.concat([all_valid, gpu_data_giulia])

# group results by: dataset, xdrop
all_devices_df = all_valid.loc[(all_valid["devices"] == 1) & (all_valid["name"].apply(lambda n: "decom" not in n or "decomn" in n))]
device_subset_df = all_devices_df.loc[all_devices_df["devicename"].isin(("AMD EPYC 7763", "Graphcore IPU", "NVIDIA A100")) & (all_devices_df["dataset"].isin(("simulated85", "celegans", "ecoli", "ecoli100")))]

# Create overview tables
overview_table = device_subset_df.loc[device_subset_df["xdrop"] <= 20][["xdrop", "dataset", "devicename", "gcups_inner", "time_ms_inner", "algo"]].set_index(["xdrop", "dataset", "devicename", "algo"]).sort_index()
overview_table_cpu = overview_table.xs("AMD EPYC 7763", level="devicename")
overview_table_ipu = overview_table.xs("Graphcore IPU", level="devicename")
overview_table_gpu = overview_table.xs("NVIDIA A100", level="devicename")

overview_table.to_csv("values.csv")

overview_ratio = (overview_table_ipu.reset_index(level="algo", drop=True) / pd.concat((overview_table_gpu, overview_table_cpu))).sort_index()
overview_ratio
overview_ratio.to_csv("ratios.csv")


# PERFORMANCE CPU - IPU Figure
fig = pplt.figure()
axes = fig.subplots(ncols=4)
artists = {}
datasets = "simulated85", "ecoli", "ecoli100", "celegans"
for i, (xdrop, xdf) in enumerate(device_subset_df.groupby("xdrop")):
    if xdrop > 20:
        continue
    ax = axes[i]
    plot_df = xdf.pivot(index="dataset", columns="algo", values="gcups_inner")
    plot_df = plot_df.reindex(datasets)
    plot_df = plot_df[["logan", "ksw2", "seqan", "ipuma"]]
    for i, ds in enumerate(datasets):
        ipuma_perf = plot_df.loc[ds, "ipuma"]
        seqan_perf = plot_df.loc[ds, "seqan"]
        offset=0.4
        ax.hlines((seqan_perf, ipuma_perf), i, i+offset - 0.2, linestyles="solid", color="red", lw=1)
        vx = i + 0.11
        ax.plot((vx, vx), (seqan_perf, ipuma_perf), c="red", ls="dashed", lw=1)
        ax.text(
            vx - 0.35, ((ipuma_perf + seqan_perf) / 2), f"{ipuma_perf / seqan_perf:.1f}×",
            color="red",
            bbox=dict(facecolor="w", ls='', pad=-.2),
        )
        # ax.axvline(vx, seqan_perf, ipuma_perf, color="red", linestyle="dashed", lw=1)

    ax.bar(plot_df, cycle="Grays", edgecolor="black")
    ax.format(ylabel="GCUPS", title=f"Xdrop = {int(xdrop)}", xlabel="", yformatter="sci")
    handles, labels = ax.get_legend_handles_labels()
    for h, l in zip(handles, labels):
        if l not in artists:
            artists[l] = h
unique = [(h, l) for i, (h, l) in enumerate(zip(handles, labels)) if l not in labels[:i]]

fig.legend(artists.values(), loc="b", ncols=4, title="Algorithm")
fig.savefig(f"datasets_xdrop_inner_total_nonlog.pdf")


# PERFORMANCE IPU by number of devices and composing/decomposing
ipu_data = all_valid.loc[all_data["device"] == "ipu"].copy()
ipu_data["multicomparison"] = ipu_data["name"].apply(lambda n: "decomn" in n)

ipu_data["dataset"].unique()
ipu_data_ds = ipu_data.loc[ipu_data["dataset"].isin(("ecoli100", "celegans"))]

# Multidevice scaling results
scaling_df = ipu_data_ds.loc[ipu_data_ds["xdrop"] == 50]
scaling_df = ipu_data_ds.sort_values(by=["devices", "multicomparison"]).copy()
scaling_df = scaling_df.set_index(["xdrop", "dataset", "devices", "multicomparison"]).sort_index()
scaling_df.columns
scaling_df = scaling_df[["time_ms", "gcups"]]
us = scaling_df.unstack()
dev1 = us.xs(1, level="devices")
for dev in (2, 4, 8, 16, 32):
    r = dev1 / us.xs(dev, level="devices")
    print(f"{dev=}")
    print(r)
# -----

# Get workload balancing results
ipu_data = all_data.loc[(all_data["device"] == "ipu") & all_data["dataset"].isin(["ecoli", "ecoli100", "celegans"]) & (all_data["devices"] == 1)].copy()
ipu_data.drop(columns=["host", "device", "devicename", "name", "algo", "threads"], inplace=True)
ipu_data.set_index(["xdrop", "devices" ,"dataset"], inplace=True)
ipu_data.sort_index(inplace=True)
indices = []
records = []
for i, s in ipu_data.iterrows():
    d = s["data"]
    joblog_data = pd.DataFrame.from_records([
        e["data"]
        for e in d["entries"]
        if e["type"] == "JOBLOG"
    ])
    indices.append(i)
    records.append({
        "transfer_info_ratio": joblog_data.mean()["transfer_info_ratio"],
        "comparison_occupancy": joblog_data.mean()["comparison_occupancy"],
        "num_batches": len(joblog_data),
    })


joblog_data[0].keys()

fig = pplt.figure()
axes = fig.subplots(ncols=5)
for i, x in enumerate((5, 10, 15, 20, 50)):
    ipu_data_x10 = ipu_data_ds.loc[ipu_data_ds["xdrop"] == x]
    ipu_data_multi = ipu_data_x10.loc[ipu_data["multicomparison"]]
    ipu_data_single = ipu_data_x10.loc[~ipu_data["multicomparison"]]
    plot_df = ipu_data_multi.pivot(index="devices", columns="dataset", values="time_ms") / 1e3
    plot_dfs = ipu_data_single.pivot(index="devices", columns="dataset", values="time_ms") / 1e3
    ax = axes[i]
    cycle1 = pplt.Cycle('blues', 2, left=0.4, right=0.8)
    cycle2 = pplt.Cycle('reds', 2, left=0.4, right=0.8)
    cycle_linestyle = {'linestyle':('-',':')}

    lines1 = ax.plot(plot_df, cycle_kw=cycle_linestyle)
    lines2 = ax.plot(plot_dfs, cycle_kw=cycle_linestyle)

    art = ax.scatter(plot_df, marker="x", linestyle="solid")
    art2 = ax.scatter(plot_dfs, marker="|", linestyle="solid")

    ax.format(yscale="log", xscale="log", ylabel="Execution time s", xlabel="Number of IPU devices", title=f"Xdrop = {x}", xlocator=pplt.Locator("fixed", ipu_data_single["devices"]))

fig.legend(
    [*lines1, art[0], art2[0]], [*[l.get_label() for l in lines1], "Multicomparison", "Singlecomparison"], loc="b", title="", ncols=4)

fig.savefig("multidevice_xdrop_e2egcups.pdf")
