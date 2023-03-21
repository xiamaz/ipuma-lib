#%%
import re
import numpy as np
import matplotlib.pyplot as plt

reg = re.compile(r"^T\[(\d+).(\d)\]: AAA: (-?\d+)$")

# to csv
def printcsv(fname):
        print("phase,tile,thread,time")
        with open(fname, "r") as fd:
                level = 0
                wasmatched = False
                for line in fd:
                        # T[1467.6]: AAA: 174
                        match = reg.match(line)
                        if not match:
                                if wasmatched:
                                        level += 1
                                        wasmatched = False
                                continue
                        wasmatched = True   
                        print(f"{level},{match.group(1)},{match.group(2)},{match.group(3)}")


# %%
def nestedlevels(fname):
        levels = []
        level = []
        wasmatched = False
        with open(fname, "r") as fd:
                for line in fd:
                        # T[1467.6]: AAA: 174
                        match = reg.match(line)
                        if not match:
                                if wasmatched:
                                        levels.append(level)
                                        level = []
                                        wasmatched = False
                                continue
                        wasmatched = True   
                        tile = int(match.group(1))
                        thread = int(match.group(2))
                        val = int(match.group(3))
                        level.append(val)
        levels.pop()
        return levels

# %%
def nestedlevels_validate(fname):
        reg = re.compile(r"^T\[(\d+).(\d)\]: AAA:")
        levels = []
        level = np.zeros((1472, 6, 20), dtype=np.int32)
        wasmatched = False
        with open(fname, "r") as fd:
                for line in fd:
                        # T[1467.6]: AAA: 174
                        match = reg.match(line)
                        if not match:
                                if wasmatched:
                                        levels.append(level)
                                        level = np.zeros((1472, 6, 20), dtype=np.int32)
                                        wasmatched = False
                                continue
                        vals = line.split("AAA:")[1].split(" ")
                        wasmatched = True   
                        tile = int(match.group(1))
                        thread = int(match.group(2)) - 1
                        for i in range(20):
                                val = int(vals[i+1])
                                level[tile, thread, i] = val
        levels.pop()
        return levels


#%%
# data = nestedlevels("/home/lukb/git/ipuma-lib/build/test/balanceX50.log")
# data = nestedlevels("/home/lukb/git/ipuma-lib/build/test/balanceX10_model1.log")
# data = nestedlevels("/home/lukb/git/ipuma-lib/build/test/balanceX10_doubleseed.log")
# data = nestedlevels("/home/lukb/git/ipuma-lib/build/test/balanceX10_singlecount.log")
# data = nestedlevels("/home/lukb/git/ipuma-lib/build/test/balanceX10_maincount.log")
# data = nestedlevels("/home/lukb/git/ipuma-lib/build/test/run.log")
data = nestedlevels("/home/lukb/git/ipuma-lib/build/test/balanceX10_perfectorder.log")
# data = nestedlevels_validate("/home/lukb/git/ipuma-lib/build/test/balanceX10_threead_occ.txt")
# data = nestedlevels("/home/lukb/git/ipuma-lib/build/test/balanceX10_newcounts.txt")
# data = nestedlevels("/home/lukb/git/ipuma-lib/build/test/balanceX10_idealcomplexity.txt")

# %%
# Figure: Large tile hist fig
fig, axs = plt.subplots(len(data)//3 + 1, 3)

for i, level in enumerate(data):
        axs[i//3, i%3].hist(level, bins=30)
        # plt.plot()
fig.set_figheight(30)
fig.set_figwidth(10)

#%%
for d in data:
        x = d.sum(axis=1)
        for tile in range(1472):
                q = x[tile]
                if not np.all(q[:-1] >= q[1:]):
                        print(f"Tile error: {tile}")
                        print(q)


#%%
x = 0
for level in data:
        s = np.sum(level)
        # print(s)
        x += s
print(x)

# %%
# Figure: Max / avg(work)
imbalances = []
plt.figure()
for level in data:
        imbalances.append(np.max(level)/(np.sum(level)/len(level)))
print(f"Average Imbalance {np.average(imbalances)}")
plt.plot(imbalances)
plt.title("Max Over Avg")
plt.ylabel("Max/avg(threads)")
plt.xlabel("Batch Number")
plt.draw()

# %%
# Figure: Utilization [%] 
utilizations = []
plt.figure()
for level in data:
        full_worktime = np.max(level) * len(level)
        used_worktime = np.sum(level)
        utilizations.append((used_worktime/full_worktime) * 100.0)
print(f"Average Utilization {np.average(utilizations)}")
plt.plot(utilizations)
plt.title("Avg WorkerCtx utilization")
plt.ylabel("Utilization [%]")
plt.xlabel("Batch Number")
plt.draw()

# =====================================SPARSE MATRIX============================================

#%%
# LEN_V = 0
# LEN_H = 1
# SEED_V = 2
# SEED_H = 3
# SCORE = 4
# TIME_LEFT = 5
# TIME_RIGHT = 6
# def parse_seedinfo(fname):
#         #AAA: 16319 17595 10463 3506 615 35661936 73393326
#         levels = []
#         level = []
#         wasmatched = False
#         with open(fname, "r") as fd:
#                 for line in fd:
#                         if not line.startswith("AAA:"):
#                                 if wasmatched:
#                                         levels.append(level)
#                                         level = []
#                                         wasmatched = False
#                                 continue
#                         wasmatched = True
#                         splits = line.split(" ")
#                         level.append(tuple([int(x) for x in splits[1:]]))
#         return levels
        
# def seedinfo_to_csv(seeddata, fname):
#         with open(fname, "w") as fd:
#                 fd.write("level,seqH,seqV,seedH,seedV,score,left,right,total\n")
#                 for i,level in enumerate(seeddata):
#                         for cmp in level:
#                                 fd.write(f"{i},{cmp[LEN_H]},{cmp[LEN_V]},{cmp[SEED_H]},{cmp[SEED_V]},{cmp[SCORE]},{cmp[TIME_LEFT]},{cmp[TIME_RIGHT]},{cmp[TIME_LEFT] + cmp[TIME_RIGHT]}\n")



# seeddata = parse_seedinfo("/home/lukb/git/ipuma-lib/build/test/comptime.log")
# seedinfo_to_csv(seeddata, "/home/lukb/git/ipuma-lib/build/test/comptime.csv")

#%%

from scipy.io import mmread

def read_values(fname):
        vtxmap = {}
        with open(fname, "r") as fd:
                for line in fd:
                        (s_vtx, s_val) = line.split(" ")
                        vtxmap[int(s_vtx)] = int(s_val)
        m = np.max(list(vtxmap.keys()))
        print(m)
        vals = np.zeros((m+1,), dtype=np.int32)
        for i in range(m+1):
                vals[i] = vtxmap.get(i, 0)
        return vals

vals = read_values("/home/lukb/git/ELBA/ELBA/build_release/val.txt")

#%%

from scipy.io import mmread
from scipy.sparse import coo_matrix

def read_coo(fname):
        row = []
        col = []
        with open(fname, "r") as fd:
                for line in fd:
                        (s_row, s_col) = line.split(" ")
                        row.append(int(s_row))
                        col.append(int(s_col))
        N = len(row)
        rowmax = np.max(row)
        colmax = np.max(col)
        data = np.ones((N,), dtype=np.int32)
        s = max(rowmax+1, colmax+1)
        print(f"Read {N} values shape=({(s, s)})")
        print(f"Avg val/[row,col] => {N/s}")
        return coo_matrix((data, (row, col)), (s, s))

A = read_coo("/home/lukb/git/ELBA/ELBA/build_release/mat2.txt")

#%%
plt.figure(figsize=(20, 20), dpi=100)
plt.spy(A.tocsc()[:2000,:2000])

#%%
import metis


# %%
# Reverse Cuthill-McKee
from scipy.sparse.csgraph import reverse_cuthill_mckee
# graph = reverse_cuthill_mckee(A.tocsc())

# Aperm = coo_matrix((A.data, (graph[A.row], graph[A.col])), A.shape)
# plt.figure(figsize=(20, 20), dpi=100)
# plt.spy(Aperm.tocsc()[:2000,:2000])

print("create sym matrix")
Ain =  A.copy().transpose() + A
print("Sym matrix created")

print("reverse_cuthill_mckee")
rcm_perm = reverse_cuthill_mckee(Ain.tocsr())
# rcm_perm = np.argsort(vals)

coo_Ain = Ain.tocoo()

print("Permute dict")
# rev_perm_dict = {k : rcm_perm.tolist().index(k) for k in rcm_perm}

rev_perm = np.zeros(rcm_perm.shape, dtype=np.int32)
for i in range(rcm_perm.shape[0]):
        rev_perm[rcm_perm[i]] = i

# y = np.arange(0, np.max(rcm_perm)+1)
# rev_perm_dict = y[rcm_perm]
print("Permute arr")
perm_i = rev_perm[coo_Ain.row]
perm_j = rev_perm[coo_Ain.col]

print("new matrix")
new_matrix = coo_matrix(
    (Ain.data, (perm_i, perm_j)), 
    shape=Ain.shape
)

#%%
plt.figure(figsize=(20, 20), dpi=100)
viewport = 200
x = 8000
y = 8000
plt.spy(new_matrix.tocsc()[y:y+viewport, x:x+viewport])

# %%

def load(fname):
        with open(fname) as fd:
                return np.array(json.load(fd)[0])


cpu = load("/tmp/parity/scores_cpu_seqan_sim05.json")
ipu_ipuint = load("/tmp/parity/scores_ipuint_sim05.json")
ipu_ipufloat = load("/tmp/parity/scores_ipu_sim05.json")
ipu_cpu = load("/tmp/parity/scores_cpu_ipumacpu_sim05.json")

print("CPU IPU")
print(st.describe(ipu_cpu/cpu))

print("CPU ipuint/ipufloat")
print(st.describe(ipu_ipufloat/ ipu_ipuint))

print("CPU ipu/ipu")
print(st.describe(ipu_ipuint/ ipu_cpu))


print("Diff blocks")

for d in [0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0]:
        print(d, np.where(np.abs((cpu/ipu_cpu) - 1) <= d)[0].size)
