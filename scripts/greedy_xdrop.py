from typing import Dict, Any
from math import inf, floor

sequences1 = [
    "AAA",
    "AAAA",
    "ATAT",
    "AAAT",
    "AAAA",
    "AAAAAAAAAAAAAAA",
    "AAATTTTTTTTTTT",
    "TTTTTTTTTTTAAA",
    # "aaaaaaaaaaaaa"
    "CATCGCCCGATTTTCACGTTCGAGAGCGGCGGAGCGGATCGCTCCTTGTTCTTTTTGCCAGGCCCGTAGTTCTTCACCCGTTTTGAATTCGGGTTTGTAT",
    "GCCAGGCAAAATCGGCGTTTCTGGCGGCGATGAGCCATGAGATCCGCACACCGCTGTACGGTATTCTCGGCACTGCTCACTTGATGGCAGATAACGCGCC",
]
sequences2 = [
    "AAA",
    "AAA",
    "ATAT",
    "AAAA",
    "AAAT",
    "AAATTTTTTTTTTTT",
    "AAA",
    "TTAAATT",
    "CATCGCCCGATTTTCACGTTCGAGAGCGGCGGAGCGGATCGCTCCTTGTTCTTTTTGCCAGGCCAGTAGTTCTTCACCCGTTTTGAATGCGGGTTTGATA",
    "GCCAGGCAAAATCGGCGTTTCTGGCGGCGATGAGCCATGAGATCCGCACACCGCTGTACGGTATTCTCGGCACTGCTCAACTGCTGGCAGATAACCCCGC",
    # "aaabbaaaaaaaa"
]
X = 10

def greedyXdrop(s1, s2):

    M = len(s1)
    N = len(s2)

    mis = -2
    mat = 2

    def calcSscore(ij, d):
        return (ij)*mat/2 - d * (mat-mis)

    # print(f"{s1=} {s2=}")
    # print(f"{M=} {N=} {calcSscore(M+N, 2)=}")

    i = 0
    while i < min(M, N) and s1[i] == s2[i]:
        i += 1

    R: Dict[tuple, Any] = {(int(0), int(0)): i}

    xdrop_offset = floor((X + mat / 2) / (mat - mis)) + 1
    T = [0.0] * (M+N + xdrop_offset+ 1)
    T[0 + xdrop_offset] = calcSscore(i + i, 0)
    Tt = T[0 + xdrop_offset]

    d = 0
    L = 0
    U = 0

    while True:
        # print("-=-----")
        d += 1
        dd = d - xdrop_offset
        for k in range(L-1, U+1 + 1):
            # print(f"{k=} {L=} {U=}")
            ivals = [-inf]
            if L < k:
                ivals.append(
                    R[(d - 1, k - 1)] + 1
                )
            if L <= k and k <= U:
                ivals.append(
                    R[(d - 1, k)] + 1
                )
            if k < U:
                ivals.append(
                    R[(d - 1, k + 1)]
                )
            # print(f"{ivals}, {i=}, {k=}, {d=} {dd=} {U=} {L=}")
            i = max(ivals)
            j = i - k
            # print(f"{ivals}, {i=}, {j=}, {k=}, {d=} {dd=}, {calcSscore(i+j, d)}")
            if i > -inf and calcSscore(i+j, d) >= T[dd + xdrop_offset] - X:
                while i < M and j < N and s1[i] == s2[j]:
                    i += 1
                    j += 1
                R[(d, k)] = i
                Tt = max(Tt, calcSscore(i+j, d))
            else:
                R[(d, k)] = -inf
            # print(R)

        # print(f"{d=} {xdrop_offset=}")
        T[d+xdrop_offset] = Tt
        L = min(k for (_, k), v in R.items() if v > -inf)
        U = max(k for (_, k), v in R.items() if v > -inf)
        # from pprint import pprint
        # print("-=-----")
        # print([k for (_, k), v in R.items() if v == M])
        # pprint({k: v for k, v in R.items() if v > -inf})
        L = max(L, max([k for (_, k), v in R.items() if v == N + k] or [-inf]) + 2)
        U = min(U, min([k for (_, k), v in R.items() if v == M] or [inf]) - 2)

        # print(f"{Tt} {L=} {U=} {M=} {N=}")

        if L > U + 1 or d >= M + N:
            break

    return Tt

for s1, s2 in zip(sequences1, sequences2):
    print(f"####\n{s1=}\n{s2=}\n   {greedyXdrop(s1, s2)=}\n####")
