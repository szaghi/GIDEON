#!/usr/bin/env python3
"""Draw the figures of GIDEON's documentation from the report data, in the light and dark themes.

Every curve and every number of a figure is read from a file: the generator's JSON (the exact rows; spectral curves
are computed from them), the CSV files and outputs of the test programs run by scripts/run_report.sh (e(h), observed
orders, the ADAM comparisons), the recorded mutation run, and the benchmark data in notes/bench-data. Each figure
<name>.png is drawn on washi and, with --theme both, again as <name>-dark.png on smoked bronze (figure_theme.py); the
docs theme swaps them.

    docs/public/guide/windows-{fd,fv,upwind}.png      the windows of each geometry, with the weights of a row
    docs/public/guide/performance.png                 GIDEON / ADAM cost, GPU and one CPU core (notes/bench-data)
    docs/public/verification/spectral-centred.png     modified wavenumber of the centred d1 and d2
    docs/public/verification/spectral-upwind.png      dispersion and dissipation of the left-side d1
    docs/public/verification/convergence.png          e(h) in R8P and R16P: the round-off floor
    docs/public/verification/orders.png               observed - design order of every kernel, R8P and R16P
    docs/public/verification/adam.png                 the comparisons with ADAM by class
    docs/public/verification/advection.png            the wave packet after two periods (gideon_vv_advection)
    docs/public/regression/mutation.png               which test kills which mutant

Usage: make_figures.py [--report DIR] [--theme light|dark|both] [--only NAME ...]
"""

from __future__ import annotations

import argparse
import csv
import re
import statistics
import sys
from collections import defaultdict
from fractions import Fraction as Fr
from pathlib import Path

import matplotlib as mpl

mpl.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap  # noqa: E402
from matplotlib.patches import Rectangle  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
import gideon_palette as P  # noqa: E402
from figure_theme import theme, themes  # noqa: E402
from make_tables import ADAM_GROUP, FAMILY, Report, short, symbol_exact  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
PUB = ROOT / 'docs' / 'public'
DPI = 150
plt.rcParams.update({'font.size': 9, 'axes.titlesize': 9, 'axes.labelsize': 9, 'legend.fontsize': 8,
                     'xtick.labelsize': 8, 'ytick.labelsize': 8, 'savefig.dpi': DPI, 'savefig.bbox': 'tight',
                     'savefig.pad_inches': 0.08, 'axes.spines.top': False, 'axes.spines.right': False})


class Ink:
    """The colours of one theme (gideon_palette)."""

    def __init__(self, dark: bool) -> None:
        self.dark = dark
        self.centred = P.CENTRED_D if dark else P.CENTRED
        self.upwind = P.UPWIND_D if dark else P.UPWIND
        self.right = P.RIGHT_D if dark else P.RIGHT
        self.adam = P.ADAM_D if dark else P.ADAM
        self.ink = 'k'  # mapped to the dark ink by figure_theme in the dark theme
        self.blue_ramp = P.RAMP_BLUE_D if dark else P.RAMP_BLUE
        self.moegi_ramp = P.RAMP_MOEGI_D if dark else P.RAMP_MOEGI

    def cmap(self, ramp: tuple[str, ...]) -> LinearSegmentedColormap:
        return LinearSegmentedColormap.from_list('ramp', list(ramp))

    def shades(self, ramp: tuple[str, ...], n: int) -> list:
        c = self.cmap(ramp)
        return [c(i / max(n - 1, 1)) for i in range(n)]


def frac(s: str) -> str:
    """A weight as a short fraction label."""
    f = Fr(s)
    return str(f.numerator) if f.denominator == 1 else f'{f.numerator}/{f.denominator}'


# ---------------------------------------------------------------------------------------------------------------------
# the windows (guide/conventions)
# ---------------------------------------------------------------------------------------------------------------------
def window_row(ax, k: dict, ink: Ink, colour, title: str) -> None:  # noqa: ANN001
    """One kernel's window on the offset axis: data (points or cells), the target, the weights of the row."""
    lo, hi = k['lo'], k['hi']
    geom = {'fd_face_c': 'fdface', 'fd_node_c': 'fdnode', 'fv_face_c': 'fvface', 'fv_facedn_c': 'fvface',
            'fv_face_u': 'fvface', 'fv_facedn_u': 'fvface', 'fv_cell_c': 'fvcell', 'fv_cell_u': 'fvcell'}[k['family']]
    w = {int(o): v for o, v in k['weights'].items()}
    ax.plot([-5.7, 6.0], [0, 0], color='0.55', lw=0.8, zorder=0)
    for x in range(-5, 6):
        ax.plot([x, x], [-0.06, 0.06], color='0.55', lw=0.8, zorder=0)
    for o in range(lo, hi + 1):
        if geom == 'fdnode':
            xc, a, b = o, None, None
        elif geom == 'fdface':
            xc, a, b = o - 0.5, None, None
        elif geom == 'fvface':
            xc, a, b = o - 0.5, o - 1, o
        else:
            xc, a, b = o, o - 0.5, o + 0.5
        if a is None:
            ax.plot([xc], [0], 'o', ms=8, color=colour, zorder=3)
        else:
            ax.add_patch(Rectangle((a + 0.04, -0.28), b - a - 0.08, 0.56, facecolor=colour, alpha=0.28,
                                   edgecolor=colour, lw=1.0, zorder=2))
            ax.plot([xc], [0], 's', ms=4, color=colour, zorder=3)
        ax.text(xc, -0.5, f'{o}', ha='center', va='top', fontsize=8, color='0.35')
        if o in w:
            ax.text(xc, 0.42, frac(w[o]), ha='center', va='bottom', fontsize=8, color=ink.ink)
    # the target
    if geom in ('fdface', 'fvface'):
        ax.plot([0, 0], [-0.42, 0.36], color=ink.ink, lw=2.0, zorder=4)
        ax.text(0, -0.92, 'face x = 0', ha='center', va='top', fontsize=7.5, color=ink.ink)
    elif geom == 'fdnode':
        ax.plot([0], [0], marker='o', ms=12, mfc='none', mec=ink.ink, mew=1.6, zorder=4)
    else:
        ax.add_patch(Rectangle((-0.5, -0.36), 1.0, 0.72, facecolor='none', edgecolor=ink.ink, lw=1.8, zorder=4))
    ax.text(-9.4, 0.0, title, ha='left', va='center', fontsize=9, color=ink.ink)
    ax.set_xlim(-9.5, 6.1)
    ax.set_ylim(-1.15, 0.9)
    ax.axis('off')


def fig_windows(r: Report, ink: Ink) -> None:
    k = r.k
    rows = {
        'windows-fd': [
            (k['fd_d1_c4'], ink.centred, '`fd_d1_c4`\nnode 0, q(-2:2)\npoints at x = o'),
            (k['fd_face_c4'], ink.centred, '`fd_face_c4`\nface 0, q(-1:2)\npoints at x = o − 1/2'),
        ],
        'windows-fv': [
            (k['fv_face_c4'], ink.centred, '`fv_face_c4`\nface 0, q(-1:2)\ncell o = [o − 1, o]'),
            (k['fv_d1_c4'], ink.centred, '`fv_d1_c4`\ncell 0, q(-2:2)\ncell o = [o − ½, o + ½]'),
        ],
        'windows-upwind': [
            (k['fv_face_ur3'], ink.right, '`fv_face_ur3`\nface 0, q(1:3)\nright side'),
            (k['fv_face_ul3'], ink.upwind, '`fv_face_ul3`\nface 0, q(-2:0)\nleft side'),
            (k['fv_d1_ur3'], ink.right, '`fv_d1_ur3`\ncell 0, q(0:3)\nright side'),
            (k['fv_d1_ul3'], ink.upwind, '`fv_d1_ul3`\ncell 0, q(-3:0)\nleft side'),
        ],
    }
    for name, spec in rows.items():
        fig, axs = plt.subplots(len(spec), 1, figsize=(8.8, 1.25 * len(spec)))
        for ax, (kk, col, title) in zip(axs, spec, strict=True):
            window_row(ax, kk, ink, col, title.replace('`', ''))
        fig.text(0.98, 0.01, 'numbers above: the weights of the row (exact, h = 1); below: the offsets',
                 ha='right', va='bottom', fontsize=7, color='0.4')
        save(fig, PUB / 'guide' / f'{name}.png')


# ---------------------------------------------------------------------------------------------------------------------
# spectral (A3)
# ---------------------------------------------------------------------------------------------------------------------
def fig_spectral(r: Report, ink: Ink) -> None:
    xi = np.linspace(1e-3, np.pi, 600)
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(8.0, 3.3))
    for ax, n in ((a1, 1), (a2, 2)):
        ks = [k for k in r.rows['kernels'] if k['family'] == 'fd_node_c' and k['n'] == n]
        cols = ink.shades(ink.blue_ramp, len(ks))
        ex = xi if n == 1 else xi**2
        ax.plot(xi / np.pi, ex, '--', color=ink.ink, lw=1.0, label='exact')
        for k, c in zip(ks, cols, strict=True):
            ax.plot(xi / np.pi, symbol_exact(k, xi).real, color=c, lw=1.3, label=f'c{k["order"]}')
        ax.set_xlabel(r'$\xi/\pi$  ($\xi = kh$)')
        ax.set_ylabel(r'Re $\tilde\xi$' if n == 1 else r'$\tilde\xi^2$')
        ax.set_title('first derivative' if n == 1 else 'second derivative')
        ax.set_xlim(0, 1)
        ax.set_ylim(0, None)
        ax.grid(alpha=0.25)
    a1.legend(loc='upper left', frameon=False)
    fig.text(0.5, -0.02, 'FD node and FV cell-average kernels of one order have the same row: one curve for both',
             ha='center', fontsize=7.5, color='0.4')
    fig.tight_layout()
    save(fig, PUB / 'verification' / 'spectral-centred.png')

    fig, (a1, a2) = plt.subplots(1, 2, figsize=(8.0, 3.3))
    ks = [k for k in r.rows['kernels'] if k['family'] == 'fv_cell_u' and k['n'] == 1 and k['bias'] == 'ul']
    cols = ink.shades(ink.moegi_ramp, len(ks))
    a1.plot(xi / np.pi, xi, '--', color=ink.ink, lw=1.0, label='exact')
    a2.axhline(0, color=ink.ink, lw=0.9, ls='--')
    for k, c in zip(ks, cols, strict=True):
        z = symbol_exact(k, xi)
        a1.plot(xi / np.pi, z.real, color=c, lw=1.2, label=f'ul{k["order"]}')
        a2.plot(xi / np.pi, z.imag, color=c, lw=1.2)
    a1.set_title('dispersion (the wide stencils leave the frame)')
    a1.set_ylim(-2, 6)
    a1.set_ylabel(r'Re $\tilde\xi$')
    a2.set_title('dissipation (Im > 0 amplifies)')
    a2.set_ylabel(r'Im $\tilde\xi$')
    a2.set_yscale('symlog', linthresh=0.01)
    for ax in (a1, a2):
        ax.set_xlabel(r'$\xi/\pi$')
        ax.set_xlim(0, 1)
        ax.grid(alpha=0.25)
    a1.legend(loc='upper left', ncol=2, frameon=False)
    fig.tight_layout()
    save(fig, PUB / 'verification' / 'spectral-upwind.png')


# ---------------------------------------------------------------------------------------------------------------------
# orders (A1)
# ---------------------------------------------------------------------------------------------------------------------
CONV = ('fd_d1_c4', 'fd_d1_c10', 'fv_face_ur9', 'fv_d6_ur4')


def read_steps(r: Report, tag: str) -> dict[str, list[dict]]:
    bits = {'r32': 32, 'r64': 64, 'r128': 128}[tag]
    out = defaultdict(list)
    for row in csv.DictReader((r.d / tag / f'order-r{bits}.csv').open()):
        out[row['kernel']].append(row)
    return out


def fig_convergence(r: Report, ink: Ink) -> None:
    steps = {tag: read_steps(r, tag) for tag in ('r64', 'r128')}
    fig, axs = plt.subplots(1, 4, figsize=(9.6, 2.9))
    styles = {'r128': (ink.centred, 'R16P'), 'r64': (ink.adam, 'R8P')}
    for ax, name in zip(axs, CONV, strict=True):
        p = int(r.k[name]['order'])
        for tag in ('r128', 'r64'):
            rows = steps[tag][name]
            h = np.array([float(x['h']) for x in rows])
            e = np.array([float(x['e']) for x in rows])
            fl = np.array([float(x['floor']) for x in rows])
            ok = np.array([x['resolved'] == 'T' for x in rows])
            col, lab = styles[tag]
            ax.loglog(h, np.maximum(e, 1e-40), '-', color=col, lw=1.0, alpha=0.7)
            ax.loglog(h[ok], e[ok], 'o', ms=4, color=col, label=f'{lab}, resolved')
            ax.loglog(h[~ok], np.maximum(e[~ok], 1e-40), 'o', ms=4, mfc='none', color=col, label=f'{lab}, not')
            ax.loglog(h, fl, ':', color=col, lw=1.0, label=f'{lab} round-off floor')
        rows = steps['r128'][name]
        s = r.summary('r128')[name]
        kb = int(s['kbest'])
        h0, e0 = float(rows[kb - 1]['h']), float(rows[kb - 1]['e'])
        hh = np.array([2.0**-1, 2.0**-20])
        ax.loglog(hh, e0 * (hh / h0) ** p, '--', color='0.5', lw=0.8, label=r'$h^{p}$, p the design order')
        ax.set_title(f'{name}, order {p}')
        ax.set_xlabel('h')
        ax.set_ylim(1e-36, None)
        ax.grid(alpha=0.2)
    axs[0].set_ylabel('|kernel − exact|')
    h, lab = axs[0].get_legend_handles_labels()
    fig.legend(h, lab, loc='lower center', ncol=4, frameon=False, bbox_to_anchor=(0.5, -0.12))
    fig.tight_layout()
    save(fig, PUB / 'verification' / 'convergence.png')


SHORT_FAMILY = {'fd_face_c': 'FD face', 'fd_node_c': 'FD node d1-d6', 'fv_face_c': 'FV face',
                'fv_facedn_c': 'FV face d1-d5', 'fv_cell_c': 'FV cell d1-d6', 'fv_face_u': 'upwind face',
                'fv_facedn_u': 'upwind face d1-d5', 'fv_cell_u': 'upwind cell d1-d6', 'deconv': 'deconvolution'}


def fig_orders(r: Report, ink: Ink) -> None:
    fams = sorted(FAMILY.items(), key=lambda x: x[1][1])
    names = [short(k['name']) for f, _ in fams for k in r.rows['kernels'] if k['family'] == f]
    fam_of = {short(k['name']): k['family'] for k in r.rows['kernels']}
    x = {n: i for i, n in enumerate(names)}
    fig, ax = plt.subplots(figsize=(9.6, 3.2))
    for tag, col, lab, mk in (('r64', ink.adam, 'R8P (measured kernels)', 'o'),
                              ('r128', ink.centred, 'R16P (every kernel)', '.')):
        s = r.summary(tag)
        xs, ys = [], []
        for n in names:
            if s[n]['status'] == 'measured':
                xs.append(x[n])
                ys.append(min(float(s[n]['p1']), float(s[n]['p2'])) - int(s[n]['design']))
        ax.plot(xs, ys, mk, ms=4 if mk == 'o' else 6, mfc='none' if mk == 'o' else col, color=col,
                label=f'{lab}: {len(xs)}')
    ax.axhline(-0.25, color=ink.ink, lw=0.9, ls='--', label='bound: design − 0.25')
    ax.axhline(0, color='0.6', lw=0.6)
    edges = [0]
    for f, _ in fams:
        edges.append(edges[-1] + sum(1 for n in names if fam_of[n] == f))
    for e in edges[1:-1]:
        ax.axvline(e - 0.5, color='0.8', lw=0.6)
    ax.set_xticks([(edges[i] + edges[i + 1]) / 2 - 0.5 for i in range(len(fams))])
    ax.set_xticklabels([SHORT_FAMILY[f] for f, _ in fams], rotation=25, ha='right')
    ax.tick_params(axis='x', length=0)
    ax.set_xlim(-1, len(names))
    ax.set_ylim(-0.3, 0.12)
    ax.set_ylabel('observed − design order')
    ax.set_title(f'the {len(names)} kernels by family: the smaller of the two orders of the finest plateau')
    ax.legend(loc='lower left', frameon=False, ncol=3)
    save(fig, PUB / 'verification' / 'orders.png')


# ---------------------------------------------------------------------------------------------------------------------
# ADAM (C5)
# ---------------------------------------------------------------------------------------------------------------------
def fig_adam(r: Report, ink: Ink) -> None:
    data = list(csv.DictReader((r.d / 'r64' / 'adam.csv').open()))
    groups = [g for _, g, _ in ADAM_GROUP]
    cnt = {g: defaultdict(int) for g in groups}
    for row in data:
        g = next(g for rx, g, _ in ADAM_GROUP if re.search(rx, row['operator']))
        cnt[g][row['class']] += 1
    classes = (('bitwise', ink.centred, 'bitwise'), ('tolerance', ink.upwind, 'within 64 eps'),
               ('defect', ink.adam, 'ADAM defect (audit A)'), ('absent', '0.65', 'ADAM returns 0 (audit C)'))
    fig, ax = plt.subplots(figsize=(8.0, 4.6))
    y = np.arange(len(groups))[::-1]
    left = np.zeros(len(groups))
    for cl, col, lab in classes:
        v = np.array([cnt[g][cl] for g in groups], float)
        ax.barh(y, v, left=left, color=col, height=0.7, label=f'{lab}: {int(v.sum())}')
        left += v
    ax.set_yticks(y)
    ax.set_yticklabels([g.replace('`', '') for g in groups], fontsize=7.5)
    ax.set_xlabel('comparisons (s = 1 … 5, every component)')
    ax.legend(loc='upper right', frameon=False)
    ax.grid(axis='x', alpha=0.25)
    save(fig, PUB / 'verification' / 'adam.png')


# ---------------------------------------------------------------------------------------------------------------------
# advection (A3)
# ---------------------------------------------------------------------------------------------------------------------
ADV = (('fv_d1_c2', 'centred c2'), ('fv_d1_c4', 'centred c4'), ('fv_d1_ul2', 'upwind ul2'))


def read_frames(path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    rows = [np.array(ln.split(), float) for ln in path.read_text().splitlines() if ln.strip()]
    x = rows[0]
    n = x.size
    data = np.array(rows[1:])
    return x, data[:, 0], data[:, 1:1 + n], data[:, 1 + n:1 + 2 * n]


def fig_advection(r: Report, ink: Ink) -> None:
    cols = (ink.right, ink.centred, ink.upwind)
    fig, ax = plt.subplots(figsize=(8.0, 2.9))
    for (name, lab), col in zip(ADV, cols, strict=True):
        x, t, u, ex = read_frames(r.d / 'r64' / f'advection-frames-{name}.dat')
        if name == ADV[0][0]:
            ax.plot(x, ex[-1], '-', color=ink.ink, lw=1.0, label='exact')
        ax.plot(x, u[-1], '-', color=col, lw=1.2, label=lab)
    ax.set_xlabel('x')
    ax.set_ylabel('u (cell averages)')
    ax.set_xlim(0, 1)
    ax.set_title(f't = {t[-1]:g}: two periods, 128 cells, the carrier at ξ = π/8')
    ax.legend(loc='upper left', ncol=4, frameon=False)
    ax.grid(alpha=0.25)
    save(fig, PUB / 'verification' / 'advection.png')


# ---------------------------------------------------------------------------------------------------------------------
# performance (guide), mutation (regression)
# ---------------------------------------------------------------------------------------------------------------------
DEVICES = (('rtx4070-nvf-oac', 'RTX 4070, nvfortran OpenACC'), ('core-gnu', 'one core, gfortran'),
           ('core-nvf', 'one core, nvfortran'))


def bench() -> tuple[list[str], dict[str, dict[str, list[float]]], dict[str, list[float]]]:
    cases: list[str] = []
    ratio = {d: defaultdict(list) for d, _ in DEVICES}
    disp = {d: [] for d, _ in DEVICES}
    for d, _ in DEVICES:
        files = sorted((ROOT / 'notes' / 'bench-data').glob(f'{d}-run-*.txt'))
        if not files:
            raise SystemExit(f'notes/bench-data: no {d}-run-*.txt')
        for f in files:
            for ln in f.read_text().splitlines():
                if ln.startswith('#') or not ln.strip():
                    continue
                if 'dispatch' in ln:
                    disp[d].append(float(ln.split()[4]))
                    continue
                m = re.match(r'(.+?)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s*$', ln)
                case = m.group(1).strip()
                if case not in cases:
                    cases.append(case)
                ratio[d][case].append(float(m.group(4)))
    return cases, ratio, disp


def fig_performance(ink: Ink) -> None:
    cases, ratio, disp = bench()
    cols = (ink.centred, ink.upwind, ink.adam)
    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(9.6, 3.4), gridspec_kw={'width_ratios': [5, 1.2]})
    x = np.arange(len(cases))
    wdt = 0.26
    for i, ((d, lab), col) in enumerate(zip(DEVICES, cols, strict=True)):
        med = [statistics.median(ratio[d][c]) for c in cases]
        lo = [med[j] - min(ratio[d][c]) for j, c in enumerate(cases)]
        hi = [max(ratio[d][c]) - med[j] for j, c in enumerate(cases)]
        ax.bar(x + (i - 1) * wdt, med, wdt, color=col, label=lab, yerr=[lo, hi], capsize=1.5,
               error_kw={'elinewidth': 0.7, 'ecolor': '0.3'})
        dm = statistics.median(disp[d])
        ax2.bar([i], [dm], 0.7, color=col, yerr=[[dm - min(disp[d])], [max(disp[d]) - dm]], capsize=1.5,
                error_kw={'elinewidth': 0.7, 'ecolor': '0.3'})
    ax.axhline(1.0, color=ink.ink, lw=0.9, ls='--')
    ax.set_xticks(x)
    ax.set_xticklabels(cases, rotation=30, ha='right')
    ax.set_ylabel('time GIDEON / time ADAM')
    ax.set_ylim(0, 1.18)
    ax.legend(loc='upper left', frameon=False, ncol=3)
    ax.grid(axis='y', alpha=0.25)
    ax.set_title('kernels against ADAM (median of four runs, range as bars)')
    ax2.axhline(1.0, color=ink.ink, lw=0.9, ls='--')
    ax2.set_xticks(range(len(DEVICES)))
    ax2.set_xticklabels(['GPU', 'gfortran', 'nvfortran'], rotation=30, ha='right')
    ax2.set_ylim(0.9, 1.25)
    ax2.set_title('fd_d1_c(6) / fd_d1_c6')
    ax2.grid(axis='y', alpha=0.25)
    fig.tight_layout()
    save(fig, PUB / 'guide' / 'performance.png')


def fig_mutation(ink: Ink) -> None:
    sys.path.insert(0, str(ROOT / 'scripts'))
    run = sorted((ROOT / 'notes').glob('mutation-run-*.txt'))[-1]
    killed = {}
    for ln in run.read_text().splitlines():
        m = re.match(r'(\S+)\s+killed by (.*)$', ln)
        if m:
            killed[m.group(1)] = [t.strip().removeprefix('gideon_test_') for t in m.group(2).split(',')]
    tests = sorted({t for v in killed.values() for t in v}, key=lambda t: -sum(t in v for v in killed.values()))
    mut = list(killed)
    fig, ax = plt.subplots(figsize=(6.4, 6.0))
    for i, m in enumerate(mut):
        for j, t in enumerate(tests):
            if t in killed[m]:
                ax.add_patch(Rectangle((j - 0.42, i - 0.42), 0.84, 0.84, color=ink.centred))
    ax.set_xlim(-0.6, len(tests) - 0.4)
    ax.set_ylim(len(mut) - 0.4, -0.6)
    ax.set_xticks(range(len(tests)))
    ax.set_xticklabels(tests, rotation=35, ha='right')
    ax.set_yticks(range(len(mut)))
    ax.set_yticklabels(mut, fontsize=7.5)
    ax.tick_params(length=0)
    for s in ax.spines.values():
        s.set_visible(False)
    ax.set_title(f'{len(mut)} mutants (rows), the tests that kill each (filled)')
    save(fig, PUB / 'regression' / 'mutation.png')


def save(fig, path: Path) -> None:  # noqa: ANN001
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path)
    plt.close(fig)


FIGURES = {
    'windows': lambda r, ink: fig_windows(r, ink),
    'spectral': lambda r, ink: fig_spectral(r, ink),
    'convergence': lambda r, ink: fig_convergence(r, ink),
    'orders': lambda r, ink: fig_orders(r, ink),
    'adam': lambda r, ink: fig_adam(r, ink),
    'advection': lambda r, ink: fig_advection(r, ink),
    'performance': lambda r, ink: fig_performance(ink),
    'mutation': lambda r, ink: fig_mutation(ink),
}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument('--report', type=Path, default=ROOT / 'exe' / 'report')
    ap.add_argument('--theme', choices=('light', 'dark', 'both'), default='both')
    ap.add_argument('--only', nargs='*')
    args = ap.parse_args()
    r = Report(args.report)
    for th in themes(args.theme):
        with theme(th):
            ink = Ink(th == 'dark')
            for name, draw in FIGURES.items():
                if args.only and name not in args.only:
                    continue
                draw(r, ink)
                print(f'  {name} ({th})')
    return 0


if __name__ == '__main__':
    sys.exit(main())
