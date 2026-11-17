# A3 Spectral properties

| ID | Group | Type | Backends | Figure | Program |
|---|---|---|---|---|---|
| A3 | Accuracy | P | the rows (exact, Python); CPU (gfortran R4P, R8P, R16P; nvfortran) for the compiled kernels | yes | `gideon_test_spectral`, `gideon_vv_advection` |

## Case and scope

An order of accuracy describes the error as $h \to 0$. At a finite resolution, what an operator does to each Fourier
mode decides: how fast it propagates a wave (dispersion) and whether it damps or amplifies it (dissipation). This case
computes the modified wavenumber of every first and second derivative at a node or over a cell: from the exact rows,
for the figures and tables, and from the compiled kernels, which must agree with them. It then advects a wave packet
with three compiled kernels, so that the two effects are seen on a solution.

## Models and numerics

| Item | Choice | Code |
|---|---|---|
| Symbol | $\sigma(\xi) = \sum_o w_o e^{\mathrm{i} o\xi}$, $\xi = kh$; for cell averages the averaging factor $\sin(\xi/2)/(\xi/2)$ of the data and of the target cancels, so node (`fd`) and cell (`fv`) operators are read alike | |
| First derivative | $\tilde\xi = -\mathrm{i}\,\sigma$, exact $\tilde\xi = \xi$. Re $\tilde\xi$ is the dispersion; for $u_t + a u_x = 0$, $a > 0$, a mode grows like $e^{a\,{\rm Im}\,\tilde\xi\,t/h}$: Im $\tilde\xi < 0$ dissipates, Im $\tilde\xi > 0$ amplifies | |
| Second derivative | $\tilde\xi^2 = -\sigma$, exact $\xi^2$ | |
| From the rows | the generator's exact weights (`--json`), in double precision | [`make_tables.py`](https://github.com/szaghi/GIDEON/blob/main/docs/scripts/make_tables.py) `symbol_exact` |
| From the kernels | each kernel on $\cos(o\xi)$ and $\sin(o\xi)$ ($h = 1$), $\sigma = k(\cos) + \mathrm{i}\,k(\sin)$, at $\xi = j\pi/256$, $j = 1 \dots 256$ | [`gideon_test_spectral.F90`](https://github.com/szaghi/GIDEON/blob/main/src/tests/gideon_test_spectral.F90) |
| Asserted (frozen) | centred kernels: $\vert{\rm Im}\vert \le 64\,\varepsilon\sum\vert c\vert$ at every $\xi$; left-side d1 of 1 and 2 cells: Im $\le$ the same; every left-side kernel the mirror of its right twin, $\sigma_{ul} = (-1)^n\,\overline{\sigma_{ur}}$. The left-side d1 that amplify a band are listed, not failed | |
| Resolution | the largest $\xi/\pi$ below which $\vert{\rm Re}\,\tilde\xi - \xi\vert/\xi < 1\%$ (the resolving efficiency of Martín et al. 2006; for d2 the same on $\tilde\xi^2$) | `resolution` |
| Advection | $u_t + u_x = 0$ on the periodic $[0, 1)$, 128 cells, cell averages, RK4 at CFL 0.25, two periods; $du/dx$ from `gideon_fv_d1_c2`, `gideon_fv_d1_c4`, `gideon_fv_d1_ul2`; the packet $e^{-((x - 1/2)/0.07)^2}\sin(16\pi x)$, its carrier at $\xi = \pi/8$; asserted: the mass kept to $10^3\varepsilon$, the energy never growing, the upwind kernel losing at least half of it, c4 more accurate than c2 | [`gideon_vv_advection.F90`](https://github.com/szaghi/GIDEON/blob/main/src/tests/gideon_vv_advection.F90) |

## Program

```bash
./exe/gideon_test_spectral
./exe/gideon_vv_advection
scripts/run_report.sh          # the CSV of the symbols and the frames of the packet, then figures and casts
```

## Results

### Centred derivatives

![Modified wavenumber of the centred first and second derivatives](/verification/spectral-centred.png)

*Re $\tilde\xi$ of the first derivatives and $\tilde\xi^2$ of the second ones, orders 2 to 10. The imaginary parts
are zero. The FD node kernels and the FV cell-average kernels of one order have the same row ([C1](./rows)): one
curve for both.*

<!--@include: ./_tables/spectral-centred.md-->

### Upwind first derivatives

![Dispersion and dissipation of the left-side first derivatives](/verification/spectral-upwind.png)

*The left-side first derivatives `gideon_fv_d1_ul1` to `_ul9` (upwind for $a > 0$): Re $\tilde\xi$ (left, clipped)
and Im $\tilde\xi$ (right, on a symmetric log scale). Above the dashed line a mode grows.*

<!--@include: ./_tables/spectral-upwind.md-->

### The compiled kernels against the rows

<!--@include: ./_tables/spectral-compiled.md-->

### A wave packet

<div class="gideon-casts">
<figure class="gideon-cast" style="--ar: 2.25">
<video class="light-only" autoplay loop muted playsinline poster="/verification/casts/advection-light.png" aria-label="Linear advection of a wave packet"><source src="/verification/casts/advection-light.webm" type="video/webm"><source src="/verification/casts/advection-light.mp4" type="video/mp4"></video>
<video class="dark-only" autoplay loop muted playsinline poster="/verification/casts/advection-dark.png" aria-label="Linear advection of a wave packet"><source src="/verification/casts/advection-dark.webm" type="video/webm"><source src="/verification/casts/advection-dark.mp4" type="video/mp4"></video>
<figcaption>The packet over two periods, every frame computed by GIDEON's compiled kernels (gideon_vv_advection, report mode).</figcaption>
</figure>
</div>

![The packet after two periods](/verification/advection.png)

*After two periods: c2 lags (its group velocity at $\xi = \pi/8$ is below 1) and leaves a train of waves behind; c4
nearly keeps up; ul2 has damped the packet and shifted it.*

<!--@include: ./_tables/spectral-advection.md-->

## Findings

- **Centred kernels do not dissipate.** Their imaginary part is zero to round-off in every kind, from the rows and from
  the compiled kernels: the pair form makes the first derivative antisymmetric and the second symmetric exactly. The
  energy of the centred runs is kept to the digits printed (RK4 alone damps it, below them). A centred scheme gets its
  dissipation, if it needs one, from elsewhere.
- **Dispersion falls with the order, but slowly.** The 1% resolution grows from order 2 to 10 (table): a wider row
  resolves more of the band, never all of it. At $\xi = \pi$ every centred first derivative returns 0 (the
  odd-even mode does not move), and the second derivatives fall short of $\pi^2$.
- **Only the one- and two-cell left-side first derivatives dissipate every mode.** From three cells on, every fully
  one-sided first derivative has a band where Im $\tilde\xi > 0$: a mode of $u_t + a u_x = 0$ in that band grows,
  whatever the time step. For `gideon_fv_d1_ul3` the band is every wavenumber below $\arccos(1/4)$ (its growth rate
  is $\tfrac{a}{h}(1-\cos\xi)^2(4\cos\xi - 1)/3$); the band of each kernel is in the table. The right-side kernels
  are their mirrors and do the same for $a < 0$. These rows are exact and verified, and they are kept as building
  blocks (the candidate stencils of an upwind-biased or WENO reconstruction); as the only spatial operator of an
  advection scheme, they are unstable.
- **The compiled kernels are the rows.** Their symbols agree with the exact rows' to the precision of the kind
  (table), and the mirrors hold in the symbols as in the values ([C3](./properties)).
- **The packet confirms both effects on a solution.** With its carrier at $\xi = \pi/8$, c2 disperses the packet into
  a trailing train and keeps its energy; ul2 keeps the mass and loses most of the energy.
