# C7 The upwind–central split

| ID | Group | Type | Backends | Figure | Program |
|---|---|---|---|---|---|
| C7 | Contracts | C | Python (exact rationals); CPU (gfortran R4P, R8P, R16P; nvfortran) | no | `scripts/gen_gideon.py`, `gideon_test_split` |

## Case and scope

A hybrid scheme moves between a central and an upwind-biased face value with a sensor. That is only clean if the
biased value **is** the central one plus a dissipation term the sensor can scale: then the central scheme is recovered
exactly at $\sigma = 0$ and the biased one at $\sigma = 1$. This case proves the identity for every biased face row in
exact arithmetic, checks the compiled kernels against it, and records where it does not hold (roadmap #2, 1.2).

## Models and numerics

| Item | Choice | Code |
|---|---|---|
| Identity | $\texttt{face\_bl}_{2k+1} = \texttt{face\_c}_{2k+2} + C_k\,\delta^{2k+1}$, $\texttt{face\_br}_{2k+1} = \texttt{face\_c}_{2k+2} - C_k\,\delta^{2k+1}$, $k = 1 \dots 4$, point values and cell averages | [`gen_gideon.py`](https://github.com/szaghi/GIDEON/blob/main/scripts/gen_gideon.py), `build()` |
| Proof | the biased row minus the centred row on the centred window, divided by the binomial $\delta^{2k+1}$ row: a single ratio, or the generator stops | |
| One-sided rows | the same tried for every `ur` row against every centred face row: `split_findings` | |
| Compiled check | 50 random windows; $c + 0 \cdot d = c$ bitwise; $\vert c + d - b\vert \le 64\,\varepsilon\,(\sum\vert b\vert + \sum\vert c\vert)\max\vert q\vert$ | [`gideon_test_split.F90`](https://github.com/szaghi/GIDEON/blob/main/src/tests/gideon_test_split.F90) |

## Program

```bash
python3 scripts/gen_gideon.py --check
./exe/gideon_test_split
```

## Results

The identities, every one proven in exact rationals:

<!--@include: ./_tables/split-identities.md-->

The one-sided rows, from the generator's findings:

<!--@include: ./_tables/split-findings.md-->

The compiled kernels:

<!--@include: ./_tables/split-builds.md-->

## Findings

- Every upwind-biased face value of order $2k+1$ is the centred value of order $2k+2$ plus one face difference, on
  point values and on cell averages: the difference of two rows exact to degree $2k$ on the same $2k+2$ points
  annihilates every polynomial of degree $\le 2k$, and the only such row is $\delta^{2k+1}$. On cell averages the
  constants are $1/12$, $-1/60$, $1/280$, $-1/1260$: the familiar split of the linear WENO face values.
- **The one-sided rows do not split** beyond one cell: `ur2` − `c4` is $[1, -9, 15, -7]/16$, not a multiple of
  $\delta^3$; only `ur1` = `c2` + $\delta^1/2$. Their windows are not symmetric about the face, so the remainder
  spans more than one difference. Roadmap #2 claimed the split for "every upwind-biased and one-sided face row":
  true for the first, false for the second.
- At $\sigma = 1$ the sum $c + d$ is a different floating-point sum from the biased kernel: equal within the round-off
  of the two, not bitwise (#2's gate asked for both ends bitwise: kept at $\sigma = 0$ only).
