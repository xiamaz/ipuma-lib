from pathlib import Path

from argparse import ArgumentParser

import pandas as pd
from .ipulog_parser import load_logfile


parser = ArgumentParser()
parser.add_argument("logfile_path", type=Path)
parser.add_argument("output_path", type=Path)

args = parser.parse_args()

logfile_path = args.logfile_path
# logfile_path = Path("./xdrop_logs/elba_ecoli.log")

logs = load_logfile(logfile_path)

job_entries = pd.DataFrame.from_records(l["data"] for l in logs if l["type"] == "JOBLOG")

job_entries.describe()

job_entries.columns

import matplotlib.pyplot as plt

title = logfile_path.stem
OUTDIR = args.output_path
OUTDIR.mkdir(exist_ok=True, parents=True)
job_entries.plot.scatter("cell_count", "time_inner", title=title)
plt.savefig(OUTDIR / "timescatter.png")
job_entries[["gcups_outer", "gcups_inner"]].plot.line(title=title, xlabel="Batch Index", ylabel="GCUPS")
plt.savefig(OUTDIR / "gcups.png")
job_entries[["cell_count"]].plot.line(title=title, xlabel="Batch Index", ylabel="Cell Count")
plt.savefig(OUTDIR / "cellcount.png")
job_entries["comparison_occupancy"]
job_entries[["comparison_occupancy", "transfer_info_ratio"]].plot.line(title=title, xlabel="Batch Index", ylabel="Ratio in Percent")
plt.savefig(OUTDIR / "transfer.png")
plt.close("all")
