#%%

from pathlib import Path
import re
import numpy as np

pattern = re.compile(r"^(?P<seq1>\w+)\W+(?P<seed1>\w+)\W+(?P<seq2>\w+)\W+(?P<seed2>\w+)\W+(?P<strand>\w+)$")
As = []
Bs = []

with open("/home/zhaom/ipuma-lib/output/logan_100k/100k.txt") as f:
       lines = f.readlines()


# Good logan
i = int(0)
for l in lines:
        p = pattern.search(l)
        g = p.groupdict()
        seq1 = g["seq1"]
        seq2 = g["seq2"]
        seed1 = int(g["seed1"])
        seed2 = int(g["seed2"])
        if g["strand"] == "n":
                s1 = seq1[seed1:seed1+17]
                s2 = seq2[seed2:seed2+17]
                As.append(g["seq1"][int(g["seed1"]):][:10])
                Bs.append(g["seq2"][int(g["seed2"]):][:10])
                As.append(g["seq1"][:int(g["seed1"])][::-1][:10])
                Bs.append(g["seq2"][:int(g["seed2"])][::-1][:10])
                print(As)
                print(Bs)
                if s1 != s2:
                    raise RuntimeError(f"{i} N s1={s1} s2={s2}")
        else:
                A = seq1
                B = "".join(map(lambda x: {"A":"T", "T":"A", "G": "C", "C": "G"}[x], g["seq2"]))[::-1]
                k = 17
                _seed2 = len(B)-seed2-k

                s1 = A[seed1:seed1+17]
                s2 = B[_seed2:_seed2+17]
                # if s1 != s2:
                #     raise RuntimeError(f"{i} !N s1={s1} s2={s2}")
                # As.append(A[seed1:])
                # Bs.append(B[_seed2:])
                # As.append(A[:seed1][::-1])
                # Bs.append(B[:_seed2][::-1])
        i += 1

# Bad logan
# c = 0
# with open("elba.ipu.test.result.ipuinputmatrix.mm", "r") as lines:
#         for l in lines:
#                 l.strip("\n")
#                 p = l.split("\t")
#                 # g = p.groupdict()
#                 if len(p) != 6:
#                         continue
#                 seq1 = p[2]
#                 seq2 = p[2+2]
#                 seed1 = int(p[2+1])
#                 seed2 = int(p[2+2+1])
# 
#                 if len(seq1) > len(seq2):
#                         seq1, seq2 = seq2, seq1
#                         seed1, seed2 = seed2, seed1 
# 
#                 if len(seq1) < 10 or len(seq2) < 10:
#                         c+=1
#                         continue
#                 if len(seq1[seed1:]) > 10 and len(seq2[seed2:]) > 10:
#                         As.append(seq1[seed1:])
#                         Bs.append(seq2[seed2:])
#                 if len(seq1[:seed1]) > 10 and len(seq2[:seed2]) > 10:
#                         As.append(seq1[:seed1][::-1])
#                         Bs.append(seq2[:seed2][::-1])
# print(f"Skipped {c} lines")
                

# la = np.asarray([len(a) for a in As])
# lb = np.asarray([len(b) for b in Bs])
# for x, y in zip(la, lb):
#         print(f"c({x:5}, {y:5}):> weight {x*y}")


#%%
# import matplotlib.pyplot as plt

# plt.hist(la * lb)
# plt.hist(lb)
# plt.hist(la)

#%%

with open("elbaAs.txt", "w") as fd:
        fd.write("\n".join(As))
with open("elbaBs.txt", "w") as fd:
        fd.write("\n".join(Bs))
