#%%

from pathlib import Path
import re

pattern = re.compile(r"^(?P<seq1>\w+)\W+(?P<seed1>\w+)\W+(?P<seq2>\w+)\W+(?P<seed2>\w+)\W+(?P<strand>\w+)$")
lines = Path("./example.txt").read_text().splitlines()
As = []
Bs = []
for l in lines:
        p = pattern.search(l)
        g = p.groupdict()
        if g["strand"] == "n":
                As.append(g["seq1"][int(g["seed1"]):])
                Bs.append(g["seq2"][int(g["seed2"]):])
                As.append(g["seq1"][:int(g["seed1"])][::-1])
                Bs.append(g["seq2"][:int(g["seed2"])][::-1])



#%%
with open("As.txt", "w") as fd:
                fd.write('\n'.join(As) + '\n')
with open("Bs.txt", "w") as fd:
                fd.write('\n'.join(Bs) + '\n')

#%%
# import matplotlib.pyplot as plt
# import numpy as np
# la = np.asarray([len(a) for a in As])
# lb = np.asarray([len(b) for b in Bs])

# for x, y in zip(la, lb):
#         print(f"c({x:5}, {y:5}):> weight {x*y}")

# # plt.hist(la * lb)
# plt.hist(la)