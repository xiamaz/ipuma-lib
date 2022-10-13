from typing import Dict, Any
from math import inf, floor
# import numpy as n

sequences1 = [
    # "AAAAAAA",
    # "AAA",
    # "AAAA",
    # "ATAT",
    # "AAAT",
    # "AAAA",
    # "AAAAAAAAAAAAAAA",
    # "AAATTTTTTTTTTT",
    # "TTTTTTTTTTTAAA",
    # # "aaaaaaaaaaaaa"
    # "CATCGCCCGATTTTCACGTTCGAGAGCGGCGGAGCGGATCGCTCCTTGTTCTTTTTGCCAGGCCCGTAGTTCTTCACCCGTTTTGAATTCGGGTTTGTAT",
    "GCCAGGCAAAATCGGCGTTTCTGGCGGCGATGAGCCATGAGATCCGCACACCGCTGTACGGTATTCTCGGCACTGCTCACTTGATGGCAGATAACGCGCC",
    "GCACCGTCCAGCCAACCGCCGAGAAGAAAAGAATGAGTGCGATACGGATCAGGATATCTACGGTTTTCTGCCCCGCGCCGTTTTGCAGCCAGTTCCAGAA",
]
sequences2 = [
    # "AAAAAATTAAAAAAA",
    # "AAA",
    # "AAA",
    # "ATAT",
    # "AAAA",
    # "AAAT",
    # "AAATTTTTTTTTTTT",
    # "AAA",
    # "TTAAATT",
    # "CATCGCCCGATTTTCACGTTCGAGAGCGGCGGAGCGGATCGCTCCTTGTTCTTTTTGCCAGGCCAGTAGTTCTTCACCCGTTTTGAATGCGGGTTTGATA",
    "GCCAGGCAAAATCGGCGTTTCTGGCGGCGATGAGCCATGAGATCCGCACACCGCTGTACGGTATTCTCGGCACTGCTCAACTGCTGGCAGATAACCCCGC",
    # "aaabbaaaaaaaa"
    "CCCCGCACCCGCAAGCCGCCGAGAAAAAAAGGATGAGGGCGATACGGATCAGGATATCTACGGTTTTCTGCCCCGCGCCGTTTTGCAGCCAGTTCCAGAA",
]
X = 10


def greedyXdrop_old(s1, s2):

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
        nL = min(k for (_, k), v in R.items() if v > -inf)
        nU = max(k for (_, k), v in R.items() if v > -inf)
        # from pprint import pprint
        # print("-=-----")
        # print([k for (_, k), v in R.items() if v == M])
        # pprint({k: v for k, v in R.items() if v > -inf})
        L = max(nL, max([k for (_, k), v in R.items() if v == N + k] or [-inf]) + 2)
        U = min(nU, min([k for (_, k), v in R.items() if v == M] or [inf]) - 2)
        # print(f"old: {L=} {U=} L:{nL}<>{max([k for (_, k), v in R.items() if v == N + k] or [-inf])}")

        # print(f"{Tt} {L=} {U=} {M=} {N=}")

        if L > U + 1 or d >= M + N:
            break

    return Tt


def greedyXdrop(s1, s2):

    M = len(s1)
    N = len(s2)

    mis = -2
    mat = 2

    def calcSscore(ij, d):
        return (ij)*mat/2 - d * (mat-mis)

    i = 0
    while i < min(M, N) and s1[i] == s2[i]:
        i += 1

    R_offset = N
    R = [[-inf]*(2*N), [-inf]*(2*N)]
    R[0][0+R_offset] = i

    xdrop_offset = floor((X + mat / 2) / (mat - mis)) + 1
    T = [0.0] * (M+N + xdrop_offset+ 1)
    T[0 + xdrop_offset] = calcSscore(i + i, 0)
    Tt = T[0 + xdrop_offset]

    d = 0
    Lmin = 0
    Lmax = -inf
    L = 0
    Umin = inf
    Umax = 0
    U = 0

    if i == N:
        Lmax = 0
    if i == M:
        Umin = 0


    while True:
        d += 1
        dd = d - xdrop_offset
        for k in range(L-1, U+1 + 1):
            ivals = [-inf]
            if L < k:
                ivals.append(
                    R[0][k-1+R_offset]+1
                )
            if L <= k and k <= U:
                ivals.append(
                    R[0][k+R_offset]+1
                )
            if k < U:
                ivals.append(
                    R[0][k+1+R_offset]
                )
            i = max(ivals)
            j = i - k
            if i > -inf and calcSscore(i+j, d) >= T[dd + xdrop_offset] - X:
                while i < M and j < N and s1[i] == s2[j]:
                    i += 1
                    j += 1
                R[1][k+R_offset] = i
                Tt = max(Tt, calcSscore(i+j, d))
            else:
                R[1][k+R_offset] = -inf
            # print(R)

        print(Tt)
        T[d+xdrop_offset] = Tt

        for rk, v in enumerate(R[1]):
            k = rk - R_offset
            if v == N + k:
                Lmax = max(Lmax, k)
            if v == M:
                Umin = min(Umin, k)
            if v > -inf:
                Umax = max(Umax, k)
                Lmin = min(Lmin, k)
        L = max(Lmin, Lmax+2)
        U = min(Umax, Umin-2)
        R = R[::-1]

        if L > U + 2 or d >= M + N:
            break

    return Tt

for s1, s2 in zip(sequences1, sequences2):
    print(f"####\n{s1=}\n{s2=}\n   {greedyXdrop_old(s1, s2)=}    {greedyXdrop(s1, s2)=}\n####")
