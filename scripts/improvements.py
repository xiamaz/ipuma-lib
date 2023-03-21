#%%
import json
import os
import re
import numpy as np

#%%
XDROPS = (
#     5,
#     10,
    15,
    20
)

INPUTS = (
#     "1",
#     "2",
#     "05",
    "15",
    "elba",
)

BASEPATH = "/home/lukb/git/ipuma-lib"

EXP = (
        "single_tile",
        "single_thread",
        # "no_dynamic_switching",
        "LR_splitting",
        "static_scheduling",
        "final_no_multi_INTcodelet",
        # "final_no_multi_nofloatmax",
        "final_no_multi_ff",
        # "final_no_multi_seqlengthcplx",
        # "final_no_multi",
)

reg = re.compile(r"JOBLOG: ({.*)")
def parselog(fname):
        joblogs = []
        with open(fname, "r") as fd:
                for line in fd:
                        # print(line)
                        m = reg.search(line)
                        if not m:
                                continue
                        job = json.loads(m.group(1))
                        joblogs.append(job)
        return joblogs
                
def logstats(joblogs):
        cycles_inner = np.zeros((len(joblogs),), dtype=np.int64)
        gcups_inner = np.zeros((len(joblogs),), dtype=np.int64)
        gcups_outer = np.zeros((len(joblogs),), dtype=np.int64)

        for i, job in enumerate(joblogs):
                cycles_inner[i] = job['cycles_inner']
                gcups_inner[i] = job['gcups_inner']
                gcups_outer[i] = job['gcups_outer']

        return (np.sum(cycles_inner),
                 np.average(gcups_inner), 
                 np.average(gcups_outer), 
                 np.mean(gcups_outer))




for input in INPUTS:
        print(f"INPUT: {input}")
        for X in XDROPS:
                print(f"X-DROP: {X}")
                lastperf = None
                baseline = None
                for ex in EXP:
                        file=f"{BASEPATH}/{ex}/{input}_{X}.log"
                        data = parselog(file)
                        st = logstats(data) 
                        if lastperf == None:
                                baseline = st[0]
                                lastperf = st[0]
                        # print(f"\t{ex:>30}: {st[0]/1e7:10.2f} {st[1]:12.2f} {st[2]:12.2f} {st[3]:12.2f} {float(lastperf)/st[0]:8.2f}x {float(baseline)/st[0]:8.2f}x")
                        print(f"\t{ex:>30} & {st[0]/1e7:10.2f} & {st[1]:12.2f} & {float(lastperf)/st[0]:8.2f}x & {float(baseline)/st[0]:8.2f}x")
                        lastperf = st[0]