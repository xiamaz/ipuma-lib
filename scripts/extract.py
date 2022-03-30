#!/usr/bin/env python3
from argparse import ArgumentParser
import pathlib
import re
import json
import numpy as np
import pandas as pd

def read_log(path):
    with path.open() as f:
        data = f.read()
    return data

parser = ArgumentParser()
parser.add_argument("dataset_dir", type=pathlib.Path)
parser.add_argument("output", type=pathlib.Path)
parser.add_argument("--move_failed", type=bool, default=False)
args = parser.parse_args()

move_failed = args.move_failed

logpaths = args.dataset_dir.glob("*.log")
logdata = {p.stem : (p, read_log(p)) for p in logpaths}


def extract_ipu_lines(logdata):
    # dataset, "ipu", deviceCount, partitioning, cells, gcups_vertex, gcups_remote, gcups_run_upload, gcups_outer_loop, dropoff, inner_max_time, inner_acc_time
    table_lines = []

    for logname, (path, logtext) in logdata.items():
        lines = logtext.split("\n")
        json_log_lines = [json.loads(l.split("IPUSWLOG")[1]) for l in lines if "IPUSWLOG" in l]
        json_logs_by_type = {
            t : [j for j in json_log_lines if j["tag"] == t] for t in {j["tag"] for j in json_log_lines}
        }
        if not json_logs_by_type.get("final_log", []):
            print("No results for", logname)
            if move_failed:
                print("Moving the file to older log.")
                new_path = path.parent / f"{path.name}.old"
                path.rename(new_path)
            continue
        

        total_work_sums_l = np.int64([np.sum(b["sum_buckets"]) for b in json_logs_by_type.get("buckets", [])])
        total_work_tops_l = np.int64([b["top"] for b in json_logs_by_type.get("buckets", [])])
        total_work_maxs_l = np.int64([np.max(b["max_buckets"]) for b in json_logs_by_type.get("buckets", [])])
        intra_tile_imbalance_l = np.float32([np.average(x) for x in [np.int64(b["min_buckets"])/np.int64(b["max_buckets"]) for b in json_logs_by_type.get("buckets", [])]])
        buckets_data = {
           "bucket_global_imbalances": total_work_maxs_l / (total_work_sums_l/1472/6),
           "bucket_local_imbalances": intra_tile_imbalance_l,
        }

        # for i in range(len(total_work_sums_l)):
        #     print(f"gloabl imbalance theory: {total_work_tops_l[i] / (total_work_sums_l[i]/1472/6)}")
        #     print(f"gloabl imbalance actual: {total_work_maxs_l[i] / (total_work_sums_l[i]/1472/6)}")
        #     print(f"local imbalance: {intra_tile_imbalance_l[i]}")
        #     print("======")

        batch_perf_logs = json_logs_by_type["batch_perf"]
        setup_log, = json_logs_by_type["run_comparison_setup"]
        final_log, = json_logs_by_type["final_log"]
        batch_keys = list(batch_perf_logs[0].keys())
        funcs = {
            "median": np.median,
            "average": np.average,
            "max": np.max,
            "min": np.min,
            "sum": np.sum
        }
        batch_stat_info = {
            f"{k}_{func}": funcs[func]([g[k] for g in batch_perf_logs])
            for k in batch_keys if k not in ["tag"]
            for func in funcs
        }
        buckets_data_stat_info = {
            f"{k}_{func}": funcs[func](buckets_data[k])
            for k in buckets_data.keys()
            for func in funcs
        }
        name_parts = list(logname.split("_"))
        dataset_info = {
            **{
                "dataset": "_".join(name_parts[2:-2]),
                "numDevices": setup_log["config"]["numDevices"],
                "numThreads": setup_log["config"]["numThreads"],
                "duplication": "dup" if setup_log["config"]["duplicateDatasets"] else "nodup",
                "remote": "remote" if setup_log["config"]["ipu"]["useRemoteBuffer"] else "stream",
                "fillAlgo": setup_log["config"]["ipu"]["fillAlgo"],
            },
            **buckets_data_stat_info,
            **batch_stat_info,
            **{f"final_{k}": v for k, v in final_log.items()}
        }
        table_lines.append(dataset_info)

    return table_lines

def extract_cpu_lines(logdata):
    table_lines = []
    for logname, (path, logtext) in logdata.items():
        cpu_threads, *dsname, numa = logname.split("_")
        dataset = "_".join(dsname)
        cpu_threads = cpu_threads.split("ssw")[1]
        numa = numa.split("numa")[1]
        json_lines = [json.loads(l.split("CPUJSONLOG")[1]) for l in logtext.split("\n") if "CPUJSONLOG" in l]
        if len(json_lines) == 0:
            print(f"{logname} failed")
            if move_failed:
                print("Moving the file to older log.")
                new_path = path.parent / f"{path.name}.old"
                path.rename(new_path)
            continue

        if numa == "0-1":
            numa = 2
        else:
            numa = 1
        json_data, = json_lines
        cells = json_data["cells"]
        gcups_inner = json_data["gcups_inner"]
        gcups_outer = json_data["gcups_outer"]
        time_inner_s = json_data["time_inner_s"]
        time_outer_s = json_data["time_outer_s"]
        line = [
                dataset, "cpu-ssw", numa, cpu_threads, cells, "", gcups_inner, gcups_inner, gcups_outer, "", "", time_inner_s
        ]
        table_lines.append(line)
    return table_lines

if "ipu" in args.dataset_dir.name:
    table_lines = extract_ipu_lines(logdata)
    table = pd.DataFrame(table_lines)
    table.to_csv(args.output)
elif "cpu" in args.dataset_dir.name:
    table_lines = extract_cpu_lines(logdata)
    for l in table_lines:
        print(";".join(map(str, l)))
else:
    raise RuntimeError(f"Unsupported dataset {args.dataset_dir}")
