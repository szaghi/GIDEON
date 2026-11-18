# Dissipation dialled by a sensor

**Problem.** A central scheme does not damp what it cannot resolve: at a discontinuity it rings, and the ringing
travels. An upwind-biased scheme damps everywhere, the smooth waves included. A hybrid wants the central scheme where
the solution is smooth and the dissipation where it is not, with a sensor between the two. GIDEON writes the
dissipation as a separate kernel, so the dial is one multiplication.

## The kernels

| Kernel | Returns | Window |
|---|---|---|
| `gideon_<data>_face_diss_c<2k>` | $\delta^{2k-1}$ at the face: the undivided odd difference, the flux of a dissipation of order $2k$ | `q(1-k:k)` |
| `gideon_<data>_face_diss_bl<p>` | the dissipation part of the left-biased face value: `face_bl<p>` − `face_c<p+1>` $= C\,\delta^{p}$ | `q(-k:k+1)`, $p = 2k+1$ |
| `gideon_<data>_face_diss_br<p>` | the same for the right-biased one, $-C\,\delta^p$ | the same |

The split is exact (proven in rationals by the generator, [C7](/verification/split)):
`gideon_fv_face_bl5` = `gideon_fv_face_c6` + `gideon_fv_face_diss_bl5`, with `gideon_fv_face_diss_bl5`
$= -\delta^5/60$ (`FV_BDISS_S2`). So

```fortran
call gideon_fv_face_c6(u(i-2:i+3), c)          ! central, order 6
call gideon_fv_face_diss_bl5(u(i-2:i+3), d)    ! the dissipation of the upwind-biased bl5 (flow to the right)
f(i) = c + sigma(i)*d                          ! sigma = 0: central (bitwise); sigma = 1: bl5; between: a blend
```

Jameson, Schmidt and Turkel's blend of the second and fourth differences is written with the two smallest
differences: $f = f_c - (\varepsilon_2\, \delta^1 - \varepsilon_4\, \delta^3)$, $\delta^1$ = `gideon_fv_face_diss_c2`,
$\delta^3$ = `gideon_fv_face_diss_c4`; the switch $\varepsilon_2$, $\varepsilon_4$ is the caller's.

## Code

A smooth bump and a square wave advected once around a periodic line ($u_t + u_x = 0$, 100 cells, SSP-RK3, CFL 0.4),
with four face values: the central `c6`; `c6` plus the whole dissipation of `bl5` (that is `bl5`); the same with the
dissipation scaled by a sensor (the normalised second difference of the averages, four times, clipped at 1); and JST
on the central `c2` (Jameson's $k_2 = 1/2$, $k_4 = 1/32$):

<<< @/examples/src/cookbook_dissipation.F90{fortran}

```text
<!--@include: ../examples/output/cookbook_dissipation.txt-->
```

- **Central.** The ringing of the square wave spreads into the whole line: the bump, far from any discontinuity, ends
  with an error 17 times the biased scheme's.
- **Biased** (`sigma = 1`). The ringing is damped, and so, slightly, is the bump.
- **Sensor.** The bump keeps the central scheme (its second differences are small, the sensor stays near 0) and is the
  most accurate of the four; the square wave gets most of the dissipation, its overshoot between the central and the
  biased scheme's. A sharper sensor, or a nonlinear value where it fires, is the next step
  ([Compact-WENO rows](./compact-weno)).
- **JST.** The second difference switched on at the square wave removes the undershoot entirely, but not the
  overshoot, the largest of the four; the bump, on a second-order central scheme with fourth-order background dissipation, is
  the least accurate. JST was built for pressure switches on shocks, not for a linear square wave: the row shows the
  blend written with GIDEON's differences, not a tuned scheme.
- Every scheme conserves the total exactly: one face value per face, a flux difference per cell.

**Reference.** Roadmap issue #2, phase 1 (1.1–1.3). The split and its constants: [Operators](/guide/operators#dissipation-and-the-upwind-central-split).
The sensor and the blend are the caller's: GIDEON provides the operators (decision D7).
